/*
 * ABCReadGeo.cpp - An Alembic geometry reader for Nuke
 *
 *  Created on: Oct 24, 2011
 *  Author: Ivan Busquets
 *
 *  Expanded by: Konstantin Hristozov
 *  Updated: Aug 17, 2012
 *
 */

#define __VER__ "2.2.0"


// std libs
#include <stdio.h>
#include <iostream>
#include <time.h>
#include <sys/time.h>

// DDImage headers
#include "DDImage/ddImageVersionNumbers.h"
#include "DDImage/SourceGeo.h"
#include "DDImage/Knobs.h"
#include "DDImage/GeoInfo.h"

#include "DDImage/TableKnobI.h"
#include "DDImage/Scene.h"


// Alembic headers
#include "Alembic/AbcCoreHDF5/All.h"
#include "Alembic/AbcGeom/All.h"
#include "Alembic/Abc/All.h"

// ABCNuke helpers
#include "ABCNuke_ArchiveHelper.h"
#include "ABCNuke_GeoHelper.h"
#include "ABCNuke_MatrixHelper.h"


#define PRINT_DEBUG

#ifdef PRINT_DEBUG
#define TRACE(txt) cout << "ABCReadGeo: " << txt << endl << flush
#else
#define TRACE(txt)
#endif

#define _FPS 24.0  // Hard code a base of 24fps. Is there a way to get this from the project settings?

using namespace DD::Image;
using namespace std;
using namespace Alembic::AbcGeom;

static const char* const interpolation_types[] = { "off", "linear" , 0};
static const char* const timing_types[] = { "original timing", "retime", 0};


static const char* nodeClass = "ABCReadGeo";
static const char* const HELP = "Alembic geometry reader";

class ABCReadGeo : public SourceGeo
{
private:
    Matrix4     _local; // local matrix that Axis_Knob fills in
    Knob*       _pAxisKnob;
    
protected:
    const char*				m_filename;
    unsigned 				m_version;
    int                 		interpolate;
    IArchive 				archive;
    Knob* 				p_tableKnob;
    Table_KnobI* 			p_tableKnobI;
    int					timing;
    int					m_first;
    int					m_last;
    float				m_frame;
    float				m_sampleFrame;
    std::vector<bool>			active_objs;
    std::vector<bool>			bbox_objs;
    bool        			m_rebuild_all;
    
    vector<Alembic::AbcGeom::IObject>            _objs;
    vector<Alembic::AbcGeom::IPolyMesh>          _meshes;
    vector< vector<Alembic::AbcGeom::IXform> >   _xforms;
    
    bool            is_topology_changing;
    bool            is_deforming;
    bool            is_attribute_changing;
    bool            poly_initialized;
    bool            file_open;
    Hash            file_hash;
    Hash            new_hash;
    struct timeval tv1, tv2;
    
    
public:
	ABCReadGeo(Node* node): SourceGeo(node) 
        {
            m_filename = "";
            m_version = 0;
            interpolate = 0;
            p_tableKnob = NULL;
            p_tableKnobI = NULL;
            _pAxisKnob = NULL;
            timing = 0;
            m_first = m_last = 1;
            m_frame = 1;
            m_sampleFrame = 1;

            m_rebuild_all = true;
            file_open = false;
            is_topology_changing = false;
            is_deforming = false;
            is_attribute_changing = false;
            poly_initialized = false;
            
            _local.makeIdentity();
            
        }

	virtual void knobs(Knob_Callback f);
	int knob_changed(DD::Image::Knob* k);
	void _validate(bool for_real);
	virtual const char* Class() const {return nodeClass;}
        virtual void init_geoinfo_parms(Scene&, GeometryList&);
        virtual void build_handles(ViewerContext* ctx);

#if kDDImageVersionInteger >= 70000   
        virtual DD::Image::Op::HandlesMode  doAnyHandles(ViewerContext* ctx) {return Op::eHandlesCooked;}
#else
        virtual bool doAnyHandles(ViewerContext* ctx) {return true;}
#endif

        static const Op::Description description;

	void updateTableKnob();
	void updateTimingKnobs();


protected:
	const char * filename () const {return m_filename;}
	virtual void create_geometry(Scene& scene, GeometryList& out);
	virtual void append(Hash& hash);
	virtual void get_geometry_hash();
        void geometry_engine(Scene& scene, GeometryList& out);


};

// *****************************************************************************
// 	VALIDATE : clamp outputContext and call _validate on base class
// *****************************************************************************

void ABCReadGeo::_validate(bool for_real)
{
    
    struct timeval tv1, tv2;
    gettimeofday(&tv1, NULL);
    
    
    // original timing
    if (knob("timing")->get_value() == 0) {
            m_sampleFrame = clamp(outputContext().frame(), m_first, m_last);
    }
    //retime
    else {
            m_sampleFrame = clamp(knob("frame")->get_value_at(outputContext().frame()), m_first, m_last) ;
    }

    SourceGeo::_validate(for_real);
    
    gettimeofday(&tv2, NULL);
    float tx =  (tv2.tv_sec + tv2.tv_usec / 1000000.0) - ( tv1.tv_sec + tv1.tv_usec / 1000000.0);
    TRACE(tx << " seconds for _validate, " << 1.0 / tx << " fps");
    
}


// *****************************************************************************
// 	KNOBS : Implement the file, timing, and table knobs
// *****************************************************************************

void ABCReadGeo::knobs(Knob_Callback f)
{
        // File knobs
	File_knob(f, &m_filename, "file", "file", Geo_File);
	Button(f, "Reload");

	// Set up the common SourceGeo knobs.
	SourceGeo::knobs(f);

	// Timing knobs
	BeginGroup(f, "TimingGroup");
	Enumeration_knob(f, &interpolate, interpolation_types, "interpolation");

	Enumeration_knob(f, &timing, timing_types, "timing");

	Int_knob(f, &m_first, "first");
	Tooltip(f, "Frame of the first sample in the Alembic archive.\n"
			"You can change this to clamp the animation to a smaller framerange\n");

	Int_knob(f, &m_last, "last");
	Tooltip(f, "Frame of the last sample in the Alembic archive.\n"
			"You can change this to clamp the animation to a smaller framerange\n");
	ClearFlags(f, Knob::STARTLINE);

	Float_knob(f, &m_frame, "frame");
	Tooltip(f, "Frame the animation will be sampled from.\n"
			"You can set this to a static frame, an animated curve,\n"
			"or an expression like 'frame/2' to effectively retime the animation\n"
			"of the Alembic archive.\n\n"
			"This value will be clamped to the 'first' and 'last' values above.");
	ClearFlags(f, Knob::SLIDER);

	EndGroup(f);

	Divider(f);

	// Object management knobs
	Button(f, "activate_selection", "Activate sel.");
	Tooltip(f, "Activate selected objects\n");

	Button(f, "deactivate_selection", "Dectivate sel.");
	Tooltip(f, "Deactivate selected objects\n");

	Button(f, "bbox_selection", "Bbox sel.");
	Tooltip(f, "Set selected objects to bbox mode\n");

	Button(f, "unbbox_selection", "UnBbox sel.");
	Tooltip(f, "Unset bbox mode for selected objects\n");

	p_tableKnob = Table_knob(f, "Obj_list", "Object list");
	SetFlags(f, Knob::STARTLINE);
	Tooltip(f, "List of objects in the Alembic archive.\n"
			"For each object, the following toggles are available:\n"
			"<b>Active:</b> Enable/disable that particular object. Disabled objects will not be read from the Alembic archive.\n"
			"<b>Bbox:</b> Choose whether to read the full geometry or just a bbox of each object.\n");

        // create new tab for transform controls
        Tab_knob(f,"Transform");

        _pAxisKnob = Axis_knob(f, &_local, "transform");
        if (_pAxisKnob != NULL)
        {
            if (GeoOp::selectable() == true)
                _pAxisKnob->enable();
            else
                _pAxisKnob->disable();
        }
        
        // about tab
        Tab_knob(f, "About");
        stringstream ss;
        ss << "Compiled on " << __DATE__ << " " << __TIME__ << "\nVersion " << __VER__ << "\nKonstantin Hristozov\nDIGITAL DOMAIN";
        Newline (f, ss.str().c_str());  

        
	// Disable/enable "frame" knob based on choice in "timing" knob
	Knob* _pTimingKnob = knob("timing");
	Knob* _pFrameKnob = knob("frame");

	if (_pTimingKnob != NULL) {
		_pFrameKnob->visible(_pTimingKnob->get_value()!=0);
	}


	if (f.makeKnobs()) {
		p_tableKnobI = p_tableKnob->tableKnob();
		p_tableKnobI->addStringColumn("name", "Obj Name", false, 216 /*column width*/);
		p_tableKnobI->addColumn("active", "Active", Table_KnobI::BoolColumn, true, 45);
		p_tableKnobI->addColumn("bbox", "BBox", Table_KnobI::BoolColumn, true, 45);

	}

	// Maintain a list of active and bbox booleans, so we don't have to use the tableKnob pointer in create_geometry
	if (p_tableKnob) {
		p_tableKnobI = p_tableKnob->tableKnob();
		int numObjs = p_tableKnobI->getRowCount();
		active_objs.resize(numObjs);
		bbox_objs.resize(numObjs);
		for (int i = 0; i < numObjs; i++) {
			active_objs[i] = p_tableKnobI->getCellBool(i,1);
			bbox_objs[i] = p_tableKnobI->getCellBool(i,2);
		}
	}

}

void ABCReadGeo::build_handles(ViewerContext* ctx)
{
    // call build_matrix_handle to multiply the context model matrix with the local matrix so the
    // nodes above it will display correctly
    build_matrix_handles(ctx, _local);
}


// *****************************************************************************
// 	KNOB_CHANGED : Callbacks for file and table knobs
// *****************************************************************************

int ABCReadGeo::knob_changed(DD::Image::Knob* k)
{
    if(k->name() == "Reload") {
        knob("file")->changed();
        return 1;
    }

    if(k->name() == "Obj_list") {
        m_rebuild_all = true;
        return 1;
    }

    if(k->name() == "activate_selection") {
        p_tableKnobI->suspendKnobChangedEvents();
        std::vector<int> selectedRows = p_tableKnobI->getSelectedRows();
        for (unsigned i = 0; i < selectedRows.size(); i++) {
                p_tableKnobI->setCellBool(selectedRows[i], 1, true);

        }
        p_tableKnobI->resumeKnobChangedEvents(true);

        return 1;
    }

    if(k->name() == "deactivate_selection") {
        p_tableKnobI->suspendKnobChangedEvents();
        std::vector<int> selectedRows = p_tableKnobI->getSelectedRows();
        for (unsigned i = 0; i < selectedRows.size(); i++) {
                p_tableKnobI->setCellBool(selectedRows[i], 1, false);

        }
        p_tableKnobI->resumeKnobChangedEvents(true);

        return 1;
    }

    if(k->name() == "bbox_selection") {
        p_tableKnobI->suspendKnobChangedEvents();
        std::vector<int> selectedRows = p_tableKnobI->getSelectedRows();
        for (unsigned i = 0; i < selectedRows.size(); i++) {
                p_tableKnobI->setCellBool(selectedRows[i], 2, true);

        }
        p_tableKnobI->resumeKnobChangedEvents(true);

        return 1;
    }

    if(k->name() == "unbbox_selection") {
        p_tableKnobI->suspendKnobChangedEvents();
        std::vector<int> selectedRows = p_tableKnobI->getSelectedRows();
        for (unsigned i = 0; i < selectedRows.size(); i++) {
                p_tableKnobI->setCellBool(selectedRows[i], 2, false);

        }
        p_tableKnobI->resumeKnobChangedEvents(true);

        return 1;
    }

    if(k->name() == "file") {
        updateTableKnob();
        updateTimingKnobs();
        m_rebuild_all = true;
        return 1;
    }
    
    if (k->name() == "selectable")
    {
        if (GeoOp::selectable() == true)
            _pAxisKnob->enable();
        else
            _pAxisKnob->disable();
        return 1;
    }
    
    return SourceGeo::knob_changed(k);
}


// ***************************************************************************************
// 	UPDATETABLEKNOB : Fill in the table knob with all geo objects from the ABC archive
// ***************************************************************************************

void ABCReadGeo::updateTableKnob()
{
        TRACE("updateTableKnob()");
	p_tableKnobI->deleteAllItems();
	p_tableKnobI->reset();

	if (filename()[0] == '\0') {
		return;
	}

	IArchive archive( Alembic::AbcCoreHDF5::ReadArchive(),
			filename(),
			Abc::ErrorHandler::kQuietNoopPolicy );

	if (!archive.valid()) {
		return;
	}

	IObject archiveTop = archive.getTop();
	std::vector<Alembic::AbcGeom::IObject>  _objs;
	getABCGeos(archiveTop, _objs);

	int obj = 0;
	for( std::vector<Alembic::AbcGeom::IObject>::const_iterator iObj( _objs.begin() ); iObj != _objs.end(); ++iObj ) {
		p_tableKnobI->addRow(obj);
		p_tableKnobI->setCellString(obj,0,iObj->getName());
		p_tableKnobI->setCellBool(obj,1,true);
		obj++;
	}
}

// *****************************************************************************
// 	UPDATETIMINGKNOBS : Fill in the frame range knobs
// *****************************************************************************

void ABCReadGeo::updateTimingKnobs()
{
        TRACE("updateTimingKnobs()");
	if (filename()[0] == '\0') {
		return;
	}

	IArchive archive( Alembic::AbcCoreHDF5::ReadArchive(),
			filename(),
			Abc::ErrorHandler::kQuietNoopPolicy );

	if (!archive.valid()) {
		return;
	}

	chrono_t firstSample = std::numeric_limits<double>::max();
	chrono_t lastSample = std::numeric_limits<double>::min();

	getABCTimeSpan(archive, firstSample, lastSample);

	knob("first")->set_value(int(firstSample *  _FPS + 0.5f));
	knob("last")->set_value(int(lastSample *  _FPS + 0.5f));

}


/*----------------------------------------------------------------------------------------------------*/

// *****************************************************************************
// 	APPEND : Build up the Op's hash
// *****************************************************************************
/* virtual */
void ABCReadGeo::append(Hash& hash)
{

	p_tableKnobI = knob("Obj_list")->tableKnob();
	Op::append(hash);

	hash.append(outputContext().frame());
	hash.append(m_filename);
	hash.append(interpolate);

	if (p_tableKnobI) {
		for (unsigned i = 0; i < active_objs.size(); i++) {
			hash.append(active_objs[i]);
			hash.append(bbox_objs[i]);
		}
	}
}

// *****************************************************************************
// 	GET_GEOMETRY_HASH : Build up the geo hashes
// *****************************************************************************
/*virtual*/
void ABCReadGeo::get_geometry_hash()
{

	SourceGeo::get_geometry_hash();

	// Group Object
	//geo_hash[Group_Object].append(m_filename);

	// Group Primitives
	geo_hash[Group_Primitives].append(m_filename);
        //geo_hash[Group_Primitives].append(m_sampleFrame);
        
	// Group Points
	//geo_hash[Group_Points].append(m_filename);
	//geo_hash[Group_Points].append(m_sampleFrame);
        geo_hash[Group_Points].append(interpolate);
	
        // Group Matrix
        geo_hash[Group_Matrix].append(m_sampleFrame);
        geo_hash[Group_Matrix].append(m_filename);

        geo_hash[Group_Matrix].append(_local.a00);
        geo_hash[Group_Matrix].append(_local.a01);
        geo_hash[Group_Matrix].append(_local.a02);
        geo_hash[Group_Matrix].append(_local.a03);
    
        geo_hash[Group_Matrix].append(_local.a10);
        geo_hash[Group_Matrix].append(_local.a11);
        geo_hash[Group_Matrix].append(_local.a12);
        geo_hash[Group_Matrix].append(_local.a13);
    
        geo_hash[Group_Matrix].append(_local.a20);
        geo_hash[Group_Matrix].append(_local.a21);
        geo_hash[Group_Matrix].append(_local.a22);
        geo_hash[Group_Matrix].append(_local.a23);
    
        geo_hash[Group_Matrix].append(_local.a30);
        geo_hash[Group_Matrix].append(_local.a31);
        geo_hash[Group_Matrix].append(_local.a32);
        geo_hash[Group_Matrix].append(_local.a33);
        
	// Group Attributes
	geo_hash[Group_Attributes].append(m_filename);
	//geo_hash[Group_Attributes].append(m_sampleFrame);

	// Hash up Table knob selections
	for (unsigned i = 0; i < active_objs.size(); i++) {
		geo_hash[Group_Primitives].append(active_objs[i]);
		geo_hash[Group_Points].append(active_objs[i]);
		geo_hash[Group_Attributes].append(active_objs[i]);
		geo_hash[Group_Primitives].append(bbox_objs[i]);
		geo_hash[Group_Points].append(bbox_objs[i]);
		geo_hash[Group_Attributes].append(bbox_objs[i]);
	}

}

/*---------------------------------------------------------------------------------------------------*/

void ABCReadGeo::init_geoinfo_parms(Scene& scene, GeometryList& out) 
{
    chrono_t curTime = m_sampleFrame  / _FPS;
    gettimeofday(&tv1, NULL);
        
    const unsigned nObjects = out.size();
    TRACE("objects: " << nObjects);
    
    for (unsigned i=0; i < nObjects; ++i) 
    {
        if (!active_objs[i])
            continue;
            
        GeoInfo& info = out[i];
        Matrix4 mat = _local * getConcatMatrix( _xforms[i], curTime, interpolate != 0 );
        //mat.makeIdentity();
            
        info.matrix      = mat; //<< Set to identity in base class implementation
        info.material    = input_iop();
#if kDDImageVersionInteger >= 63100
        info.materialContext = outputContext(); // << only for 6.3+
#endif
        //
        info.render_mode = (RenderMode) render_mode_;
        info.display3d   = (Display3DMode) display3d_;
        info.selectable  = selectable();
        info.selected    = false;
        //
        info.source_geo  = 0;
        info.select_geo  = (selectable())?this:0;
        info.recursion_geo = 0;
    }
    
    gettimeofday(&tv2, NULL);
    float tx =  (tv2.tv_sec + tv2.tv_usec / 1000000.0) - ( tv1.tv_sec + tv1.tv_usec / 1000000.0);
    TRACE(tx << " seconds for getConcatMatrix, " << nObjects << " objects");
        
}
void ABCReadGeo::geometry_engine(Scene& scene, GeometryList& out)
{
    struct timeval tv1, tv2;
    gettimeofday(&tv1, NULL);
    
    // multiply the node matrix
    //for (unsigned i=0; i < out.size(); i++)
    //    out[i].matrix = _local * out[i].matrix;
    
    SourceGeo::geometry_engine(scene, out);
    
    gettimeofday(&tv2, NULL);
    float tx =  (tv2.tv_sec + tv2.tv_usec / 1000000.0) - ( tv1.tv_sec + tv1.tv_usec / 1000000.0);
    TRACE(tx << " seconds for geometry_engine, " << 1.0 / tx << " fps");
        
}

// *****************************************************************************
// 	CREATE_GEOMETRY : The meat. Query the ABC archive for the needed bits
// *****************************************************************************
/*virtual*/
void ABCReadGeo::create_geometry(Scene& scene, GeometryList& out)
{

        struct timeval tv1, tv2;
        gettimeofday(&tv1, NULL);

	if (filename()[0] == '\0') {
		out.delete_objects();
		return;
	}



        /*
	IArchive archive( Alembic::AbcCoreHDF5::ReadArchive(),
			filename(),//archiveName,
			Abc::ErrorHandler::kQuietNoopPolicy );

	if (!archive.valid()) {
		std::cout << "error reading archive" << std::endl;
		error("Unable to read file");
		return;
	}

	IObject archiveTop = archive.getTop();

	std::vector<Alembic::AbcGeom::IObject>  _objs;

	getABCGeos(archiveTop, _objs);
*/
        new_hash.reset();
        new_hash.append(filename());
        if (new_hash != file_hash || m_rebuild_all)
        {
            file_hash = new_hash;
            file_open = false;
            
            is_topology_changing = false;
            is_deforming = false;
            is_attribute_changing = false;
            
            poly_initialized = false;
            
            _objs.clear();
            _xforms.clear();
            _meshes.clear();
        }
        
        if (!file_open)
        {
            TRACE("reading file");
            
            Alembic::Abc::IArchive archive(Alembic::AbcCoreHDF5::ReadArchive(),
                                           filename(), Alembic::Abc::ErrorHandler::Policy(),
                                           Alembic::AbcCoreAbstract::ReadArraySampleCachePtr());
            
            /*
            IArchive archive( Alembic::AbcCoreHDF5::ReadArchive(),
                            filename(),//archiveName,
                            Abc::ErrorHandler::kQuietNoopPolicy );
            */
            if (!archive.valid()) 
            {
                std::cout << "error reading archive" << std::endl;
                error("Unable to read file");
                return;
            }
    
            IObject archiveTop = archive.getTop();
    
    
            //getABCGeos(archiveTop, _objs);
            getABCGeosXforms(archiveTop, _objs, _meshes, _xforms);
            
            TRACE("objs: " << _objs.size() << ", meshes: " << _meshes.size() << ", xforms: " << _xforms.size());
            
            // redo primitives if topology is changing
            for (unsigned i=0; i<_meshes.size(); i++) {
                is_topology_changing |= isTopologyChanging(_meshes[i]);
                is_deforming |= isDeforming(_meshes[i]);
                is_attribute_changing |= isNormalChanging(_meshes[i]) || isUVChanging(_meshes[i]);
            }
            
            if (is_topology_changing)
                TRACE("WARNING: topology is changing");
            if (is_deforming)
                TRACE("WARNING: deforming geo");

            
            file_open = true;
        }
        
                    

	// current Time to sample from
	chrono_t curTime = m_sampleFrame  / _FPS;


        if ( is_topology_changing )
            set_rebuild(Mask_Primitives); 
        
        if ( is_deforming )
            set_rebuild(Mask_Points);
         
        if (is_attribute_changing)
            set_rebuild(Mask_Attributes);
            
        
            
        if ( rebuild(Mask_Primitives) ) {
           TRACE("delete all objects");
	   out.delete_objects();
	}

	int obj = 0;
        for( unsigned i=0; i < _meshes.size(); i++) 
        {
            // Leave an empty obj if knob is unchecked
            if ( rebuild(Mask_Primitives) || is_topology_changing ) {
                if (!active_objs[obj] ) {
                    out.add_object(obj);
                    PointList& points = *out.writable_points(obj);
                    points.resize(0);
                    out[obj].delete_group_attribute(Group_Vertices,kUVAttrName, VECTOR4_ATTRIB);
                    obj++;
                    continue;
                }
            }

            
            // rebuild primitives if topology is changing
            if ( rebuild(Mask_Primitives) ) 
            {
                TRACE("    for loop : rebuild(Mask_Primitives): " << i);

                out.add_object(obj);
                /*
                // TODO: bbox not working now. fix it.
                if (bbox_objs[obj]) { //(bbox_mode) {
                    buildBboxPrimitives(out, obj);
                }
                else*/ 
                {
                    //buildABCPrimitives(out, obj, *iObj, curTime);
                    buildABCPrimitives(out, obj, _meshes[i], curTime);
                    
                    if (!is_deforming)
                    {
                        // set initial point positions
                        PointList& points = *out.writable_points(obj);
                        readPoly(_meshes[i], points, curTime, false);
                    }

                    poly_initialized = true;
                }
                set_rebuild(Mask_Points);
                set_rebuild(Mask_Attributes);
                TRACE("    for loop : rebuild(Mask_Primitives) END");
            }



            if ( rebuild(Mask_Points) ) 
            {
                TRACE("    for loop : rebuild(Mask_Points)");
        
                PointList& points = *out.writable_points(obj);
                /*
                // TODO: bbox not working now. fix it.

                if (bbox_objs[obj]) { 
                    Imath::Box3d bbox = getBounds(_objs[i], curTime);
    
                    points.resize(8);
                    
                    // Add bbox corners
                    for (unsigned i = 0; i < 8; i++) {
                        Vector3 pt((i&4)>>2 ? bbox.max.x : bbox.min.x, (i&2)>>1 ? bbox.max.y : bbox.min.y, (i%2) ? bbox.max.z : bbox.min.z );
                        //points[i] = xf.transform(pt);
                        points[i] = pt;
                    }
                }
        
                else*/
                {
                    struct timeval tv1, tv2;
                    gettimeofday(&tv1, NULL);
                    
                    //writePoints(*iObj, points, curTime, interpolate !=0);
                    
                    //writePoints(_meshes[i], points, curTime, interpolate !=0);
                    
                    readPoly(_meshes[i], points, curTime, poly_initialized);
                    
                    gettimeofday(&tv2, NULL);
                    float tx =  (tv2.tv_sec + tv2.tv_usec / 1000000.0) - ( tv1.tv_sec + tv1.tv_usec / 1000000.0);
                    TRACE(tx << " seconds for writePoints()");
                }
                TRACE("    for loop : rebuild(Mask_Points) END");
            }
            


		if ( rebuild(Mask_Attributes)) {
                        TRACE("    for loop : rebuild(Mask_Attributes)");
                        // TODO: bbox not working now. fix it.
                        //if (bbox_objs[obj]) { //(bbox_mode)
			//	out[obj].delete_group_attribute(Group_Vertices,kUVAttrName, VECTOR4_ATTRIB);
			//}
			//else 
                        {
				// set UVs
				Attribute* UV = out.writable_attribute(obj, Group_Vertices, kUVAttrName, VECTOR4_ATTRIB);
                                IV2fGeomParam uvParam = getUVsParam(_meshes[i]);
                                setUVs(out[obj], uvParam, UV, curTime);

				// set Normals
                                IN3fGeomParam nParam = getNsParam(_meshes[i]);
                                
        			Attribute* N = out.writable_attribute(obj, Group_Vertices, kNormalAttrName, NORMAL_ATTRIB);
				setNormals(out[obj], nParam, N, curTime);
                        }
                TRACE("    for loop : rebuild(Mask_Attributes) END");
                }
		obj++;
	}
        
	m_rebuild_all = false;
        
	out.synchronize_objects();
        
        gettimeofday(&tv2, NULL);
        float tx =  (tv2.tv_sec + tv2.tv_usec / 1000000.0) - ( tv1.tv_sec + tv1.tv_usec / 1000000.0);
        TRACE(tx << " seconds for create_geometry, " << 1.0 / tx << " fps");       
}


static Op* Build(Node* node) { return new ABCReadGeo(node); }
const Op::Description ABCReadGeo::description(nodeClass, Build);


