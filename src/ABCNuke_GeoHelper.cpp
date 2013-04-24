/*
Copyright (c) 2011, Ivan Busquets
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
 * Neither the name of Ivan Busquets nor the names of its contributors may
   be used to endorse or promote products derived from this software without
   specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
*/

//-*****************************************************************************
#include "ABCNuke_GeoHelper.h"
#include "ABCNuke_MatrixHelper.h"
#include "ABCNuke_Interpolation.h"
#include "DDImage/GeometryList.h"

// Alembic headers
#include "Alembic/AbcCoreHDF5/All.h"
#include "Alembic/AbcGeom/All.h"
//-*****************************************************************************

using namespace DD::Image;
using namespace Alembic::AbcGeom;

void fillPoints(PointList& oPointArray,
                    Alembic::Abc::P3fArraySamplePtr iPoints,
                    Alembic::Abc::P3fArraySamplePtr iCeilPoints, double alpha)
{
    //printf("fillPoints()\n");
    
    unsigned int numPoints = static_cast<unsigned int>(iPoints->size());
    oPointArray.resize(numPoints);

    if (alpha == 0 || iCeilPoints == NULL)
    {
        for (unsigned int i = 0; i < numPoints; ++i)
        {
            oPointArray[i].set( (*iPoints)[i].x, 
                                (*iPoints)[i].y, 
                                (*iPoints)[i].z);
        }
    }
    else
    {
        for (unsigned int i = 0; i < numPoints; ++i)
        {
            oPointArray[i].set( simpleLerp<float>(alpha,
                                        (*iPoints)[i].x, (*iCeilPoints)[i].x),
                                simpleLerp<float>(alpha,
                                        (*iPoints)[i].y, (*iCeilPoints)[i].y),
                                simpleLerp<float>(alpha,
                                        (*iPoints)[i].z, (*iCeilPoints)[i].z));
        }
    }

}

void readPoly(Alembic::AbcGeom::IPolyMesh iPoly, PointList& oPoints, chrono_t iFrame, bool poly_initialized)
{
    Alembic::AbcGeom::IPolyMeshSchema schema = iPoly.getSchema();
    Alembic::AbcGeom::MeshTopologyVariance ttype = schema.getTopologyVariance();
    
    Alembic::AbcCoreAbstract::index_t index, ceilIndex;
    double alpha = getWeightAndIndex(iFrame, schema.getTimeSampling(), schema.getNumSamples(), index, ceilIndex);

    //vector<Vector3> pointArray;
    Alembic::Abc::P3fArraySamplePtr ceilPoints;
    
    //printf(" * mesh samples: %d\n", (unsigned) schema.getNumSamples());
    //printf(" * normal samples: %d\n", (unsigned) schema.getNormalsParam().getNumSamples());
    //printf(" * UV samples: %d\n", (unsigned) schema.getUVsParam().getNumSamples());

    if (poly_initialized && schema.getNumSamples() <= 1)
        return;
        
    // we can just read the points
    if (ttype != Alembic::AbcGeom::kHeterogenousTopology && poly_initialized && schema.getNumSamples() > 1)
    {
        
        // read point positions into array
        Alembic::Abc::P3fArraySamplePtr points = schema.getPositionsProperty().getValue(Alembic::Abc::ISampleSelector(index));

        // if sub-frame sample, read upper frame points too
        if (alpha != 0.0)
            ceilPoints = schema.getPositionsProperty().getValue( Alembic::Abc::ISampleSelector(ceilIndex) );
        //else
        //    return;
                    
        fillPoints(oPoints, points, ceilPoints, alpha);  // copy point positions into array (interpolate if alpha != 0)
        /*
        ioMesh.setPoints(pointArray, MSpace::kObject);      // set mesh point positions
        setColors(iFrame, ioMesh, iNode.mC3s, iNode.mC4s, !iInitialized); // set colors

        // update normals, if more than one sample
        if (schema.getNormalsParam().getNumSamples() > 1)
        {
            setPolyNormals(iFrame, ioMesh, schema.getNormalsParam());
        }

        // update uvs, if more than one sample
        if (schema.getUVsParam().getNumSamples() > 1)
        {
            setUVs(iFrame, ioMesh, schema.getUVsParam());
        }
        */
    }
    else
    {
        // we need to read the topology
        Alembic::AbcGeom::IPolyMeshSchema::Sample samp;
        schema.get(samp, Alembic::Abc::ISampleSelector(index));
    
        if (alpha != 0.0 && ttype != Alembic::AbcGeom::kHeterogenousTopology)
        {
            ceilPoints = schema.getPositionsProperty().getValue(
                    Alembic::Abc::ISampleSelector(ceilIndex) );
        }
    
        fillPoints(oPoints, samp.getPositions(), ceilPoints, alpha);
    
        /*
        fillTopology(ioMesh, iParent, pointArray, samp.getFaceIndices(),
                     samp.getFaceCounts());
    
        setPolyNormals(iFrame, ioMesh, schema.getNormalsParam());
        setUVs(iFrame, ioMesh, schema.getUVsParam());
        setColors(iFrame, ioMesh, iNode.mC3s, iNode.mC4s, !iInitialized);
        */
    }
}


//////////////////////////////////////////////////////////////////////////////////////
            
void writePoints(Alembic::AbcGeom::IPolyMesh iPoly, PointList& points, chrono_t curTime = 0, bool interpolate = false) {

	IPolyMeshSchema mesh = iPoly.getSchema();
	TimeSamplingPtr ts = mesh.getTimeSampling();

	IPolyMeshSchema::Sample mesh_samp;
	const ISampleSelector iss(curTime);
	mesh.get( mesh_samp, iss );

	//IObject iObj = mesh.getObject();
	//Matrix4 xform = getConcatMatrix( iObj, curTime, interpolate );

	Vector3 pos;

        
        unsigned numPoints =  mesh_samp.getPositions()->size();
        
        
	points.resize(numPoints);
        

	Alembic::AbcCoreAbstract::index_t floorIdx = 0;
	Alembic::AbcCoreAbstract::index_t ceilIdx = 0;
	double amt = 0;

	if (interpolate) {  // check if interpolation is really needed

		amt = getWeightAndIndex(curTime, ts,
				mesh.getNumSamples(), floorIdx, ceilIdx);

		if (amt == 0 || floorIdx == ceilIdx) {
			interpolate = false;
		}
	}

	if (interpolate) {  // now get samples for interpolation

		const ISampleSelector iss_start(floorIdx);
		const ISampleSelector iss_end(ceilIdx);

		IPolyMeshSchema::Sample mesh_samp_start;
		IPolyMeshSchema::Sample mesh_samp_end;
		mesh.get( mesh_samp_start, iss_start );
		mesh.get( mesh_samp_end, iss_end );


		for (unsigned i = 0; i < numPoints; i++) {

			Imath::V3f p_start = mesh_samp_start.getPositions()->get()[i];
			Imath::V3f p_end = mesh_samp_end.getPositions()->get()[i];
			pos = lerp(Vector3(p_start.getValue()) , Vector3(p_end.getValue()), amt );

                        //points[i] = xform.transform(pos);
                        points[i] = pos;
                }

	}

	else { //no interpolation needed

		// Add points
		for (unsigned i = 0; i < numPoints; i++) {
			Imath::V3f p = mesh_samp.getPositions()->get()[i];
			pos = Vector3(p.x, p.y, p.z);
                        //points[i] = xform.transform(pos);
                        points[i] = pos;
                }
	}

}

//-*****************************************************************************

void writePoints(Alembic::AbcGeom::ISubD iSub, PointList& points, chrono_t curTime = 0, bool interpolate = false) {


	ISubDSchema mesh = iSub.getSchema();
	TimeSamplingPtr ts = mesh.getTimeSampling();
	ISubDSchema::Sample mesh_samp;

	const ISampleSelector iss(curTime);

	mesh.get( mesh_samp, iss );

	IObject iObj = mesh.getObject();
	Matrix4 xform = getConcatMatrix( iObj, curTime, interpolate );

	Vector3 pos;

	unsigned numPoints =  mesh_samp.getPositions()->size();
	points.resize(numPoints);



	Alembic::AbcCoreAbstract::index_t floorIdx = 0;
	Alembic::AbcCoreAbstract::index_t ceilIdx = 0;
	double amt = 0;

	if (interpolate) {  // check if interpolation is really needed

		amt = getWeightAndIndex(curTime, ts,
				mesh.getNumSamples(), floorIdx, ceilIdx);

		if (amt == 0 || floorIdx == ceilIdx) {
			interpolate = false;
		}
	}

	if (interpolate) {  // now get samples for interpolation

		const ISampleSelector iss_start(floorIdx);
		const ISampleSelector iss_end(ceilIdx);

		ISubDSchema::Sample mesh_samp_start;
		ISubDSchema::Sample mesh_samp_end;
		mesh.get( mesh_samp_start, iss_start );
		mesh.get( mesh_samp_end, iss_end );


		for (unsigned i = 0; i < numPoints; i++) {

			Imath::V3f p_start = mesh_samp_start.getPositions()->get()[i];
			Imath::V3f p_end = mesh_samp_end.getPositions()->get()[i];
			pos = lerp(Vector3(p_start.getValue()) , Vector3(p_end.getValue()), amt );

			points[i] = xform.transform(pos);
		}

	}

	else { //no interpolation needed

		// Add points
		for (unsigned i = 0; i < numPoints; i++) {
			Imath::V3f p = mesh_samp.getPositions()->get()[i];
			pos = Vector3(p.x, p.y, p.z);
			points[i] = xform.transform(pos);
		}
	}

}

//-*****************************************************************************

void writePoints(const Alembic::AbcGeom::IObject iObj, PointList& points, chrono_t curTime = 0, bool interpolate = false) {

	if (Alembic::AbcGeom::IPolyMesh::matches(iObj.getHeader())) {

		// Do PolyMesh
		IPolyMesh iPoly(iObj, Alembic::Abc::kWrapExisting);
		writePoints(iPoly, points, curTime, interpolate);
	}

	else if (Alembic::AbcGeom::ISubD::matches(iObj.getHeader())) {

		// Do SubD
		ISubD iSub(iObj, Alembic::Abc::kWrapExisting);
		writePoints(iSub, points, curTime, interpolate);
	}
}

//-*****************************************************************************

void fillPrimitiveIndices(const Alembic::AbcGeom::IObject iObj,
		Int32ArraySamplePtr& _fc,
		Int32ArraySamplePtr& _fi,
		chrono_t curTime = 0)
{
	if (Alembic::AbcGeom::IPolyMesh::matches(iObj.getHeader())) {
		IPolyMesh iPoly(iObj, Alembic::Abc::kWrapExisting);
		IPolyMeshSchema mesh = iPoly.getSchema();
		TimeSamplingPtr ts = mesh.getTimeSampling();
		IPolyMeshSchema::Sample mesh_samp;
		const ISampleSelector iss(curTime);
		mesh.get( mesh_samp, iss );
		_fc = mesh_samp.getFaceCounts();
		_fi = mesh_samp.getFaceIndices();
	}

	else if (Alembic::AbcGeom::ISubD::matches(iObj.getHeader())) {
		ISubD iSub(iObj, Alembic::Abc::kWrapExisting);
		ISubDSchema mesh = iSub.getSchema();
		TimeSamplingPtr ts = mesh.getTimeSampling();
		ISubDSchema::Sample mesh_samp;
		const ISampleSelector iss(curTime);

		mesh.get( mesh_samp, iss );
		_fc = mesh_samp.getFaceCounts();
		_fi = mesh_samp.getFaceIndices();
	}
}

void fillPrimitiveIndices(const Alembic::AbcGeom::IPolyMesh poly,
                          Int32ArraySamplePtr& _fc,
                          Int32ArraySamplePtr& _fi,
                          chrono_t curTime = 0)
{
    IPolyMeshSchema mesh = poly.getSchema();
    TimeSamplingPtr ts = mesh.getTimeSampling();
    IPolyMeshSchema::Sample mesh_samp;
    const ISampleSelector iss(curTime);
    mesh.get( mesh_samp, iss );
    _fc = mesh_samp.getFaceCounts();
    _fi = mesh_samp.getFaceIndices();
}

//-*****************************************************************************

Alembic::AbcGeom::IV2fGeomParam getUVsParam(const Alembic::AbcGeom::IObject iObj)
{

	IV2fGeomParam uvParam;
	if (Alembic::AbcGeom::IPolyMesh::matches(iObj.getHeader())) {
		IPolyMesh iPoly(iObj, Alembic::Abc::kWrapExisting);
		IPolyMeshSchema mesh = iPoly.getSchema();
		uvParam = mesh.getUVsParam();
	}
	if (Alembic::AbcGeom::ISubD::matches(iObj.getHeader())) {
		ISubD iSubD(iObj, Alembic::Abc::kWrapExisting);
		ISubDSchema mesh = iSubD.getSchema();
		uvParam = mesh.getUVsParam();
	}
	return uvParam;

}

Alembic::AbcGeom::IV2fGeomParam getUVsParam(const Alembic::AbcGeom::IPolyMesh poly)
{
    IV2fGeomParam uvParam;
    IPolyMeshSchema mesh = poly.getSchema();
    uvParam = mesh.getUVsParam();
    return uvParam;
}

//-*****************************************************************************

Alembic::AbcGeom::IN3fGeomParam getNsParam(const Alembic::AbcGeom::IObject iObj)
{

	IN3fGeomParam nParam;
	if (Alembic::AbcGeom::IPolyMesh::matches(iObj.getHeader())) {
		IPolyMesh iPoly(iObj, Alembic::Abc::kWrapExisting);
		IPolyMeshSchema mesh = iPoly.getSchema();
		nParam = mesh.getNormalsParam();
	}

	return nParam;

}

Alembic::AbcGeom::IN3fGeomParam getNsParam(const Alembic::AbcGeom::IPolyMesh poly)
{
    IN3fGeomParam nParam;
    IPolyMeshSchema mesh = poly.getSchema();
    nParam = mesh.getNormalsParam();

    return nParam;
}

//-*****************************************************************************

void setUVs(GeoInfo& obj,
		Alembic::AbcGeom::IV2fGeomParam & iUVs,
		Attribute* UV,
		chrono_t curTime = 0)
{

	if (!iUVs.valid())
		return;

	unsigned int numFaceVertices = 0;
	unsigned int numPrimitives = obj.primitives();
	for (unsigned p = 0; p < numPrimitives; p++) {
		for (unsigned v = 0; v < obj.primitive(p)->vertices(); v++) {
			numFaceVertices++;
		}
	}

	unsigned int numPoints = obj.points();

	Alembic::AbcGeom::IV2fGeomParam::Sample samp = iUVs.getIndexedValue();

	Alembic::AbcGeom::V2fArraySamplePtr uvPtr = samp.getVals();
	Alembic::Abc::UInt32ArraySamplePtr indexPtr = samp.getIndices();

	if (numFaceVertices != indexPtr->size() &&  numPoints != indexPtr->size()) { // UVs size is not per-point or per vertex-per-face
		return;
	}


	int uvIndex = 0;

	// per-primitive per-vertex
	if (numFaceVertices == indexPtr->size())
	{
		UV->resize(numFaceVertices);
		Vector4 _uv(0,0,0,1);

		for (unsigned pIndex = 0; pIndex < numPrimitives; ++pIndex)
		{
			int numPrimitiveVertices = obj.primitive(pIndex)->vertices();

			if (numPrimitiveVertices == 0)
				continue;

			int startPoint = uvIndex + numPrimitiveVertices - 1;  // to match the reversed winding order

			for (int vertexIndex = 0; vertexIndex < numPrimitiveVertices; vertexIndex++)
			{
				V2f uv2 = (*uvPtr)[(*indexPtr)[startPoint - vertexIndex]];
				_uv.x = uv2[0];
				_uv.y = uv2[1];
				UV->vector4(uvIndex++) = _uv;
			}
		}
	}

}

//-*****************************************************************************

void setNormals(GeoInfo& obj,
		Alembic::AbcGeom::IN3fGeomParam & Ns,
		Attribute* N,
		chrono_t curTime = 0)
{
	if (!Ns.valid())
		return;

	unsigned int numFaceVertices = 0;
	unsigned int numPrimitives = obj.primitives();
	unsigned int numPoints = obj.points();

	for (unsigned p = 0; p < numPrimitives; p++) {
		for (unsigned v = 0; v < obj.primitive(p)->vertices(); v++) {
			numFaceVertices++;
		}
	}

	Alembic::AbcGeom::IN3fGeomParam::Sample samp = Ns.getIndexedValue();

	Alembic::AbcGeom::N3fArraySamplePtr nPtr = samp.getVals();
	Alembic::Abc::UInt32ArraySamplePtr indexPtr = samp.getIndices();

	if (numFaceVertices != indexPtr->size() &&  numPoints != indexPtr->size()) { // UVs size is not per-point or per vertex-per-face
		return;
	}


	int nIndex = 0;


	// per-primitive per-vertex
	if (numFaceVertices == indexPtr->size())
	{

		N->resize(numFaceVertices);
		for (unsigned pIndex = 0; pIndex < numPrimitives; ++pIndex)
		{
			int numPrimitiveVertices = obj.primitive(pIndex)->vertices();

			if (numPrimitiveVertices == 0)
				continue;

			int startPoint = nIndex + numPrimitiveVertices - 1;  // to match the reversed winding order

			for (int vertexIndex = 0; vertexIndex < numPrimitiveVertices; vertexIndex++)
			{

				V3f normal = (*nPtr)[(*indexPtr)[startPoint - vertexIndex]];
				N->normal(nIndex++) = Vector3(normal.x, normal.y, normal.z);
			}
		}



	}


}


//-*****************************************************************************


Box3d getBounds( IObject iObj, chrono_t curTime = 0 )
{
	Box3d bnds;
	bnds.makeEmpty();

	const ISampleSelector iss(curTime);


	if ( IPolyMesh::matches( iObj.getMetaData() ) )
	{
		IPolyMesh mesh( iObj, kWrapExisting );
		IPolyMeshSchema ms = mesh.getSchema();

		bnds = ms.getSelfBoundsProperty().getValue(iss);

	}
	else if ( ISubD::matches( iObj.getMetaData() ) )
	{
		ISubD mesh( iObj, kWrapExisting );
		ISubDSchema ms = mesh.getSchema();
		bnds = ms.getSelfBoundsProperty().getValue(iss);
	}


	return bnds;
}

//-*****************************************************************************

bool isTopologyChanging(IObject iObj)
{
    if (Alembic::AbcGeom::IPolyMesh::matches(iObj.getHeader())) {
        IPolyMesh iPoly(iObj, Alembic::Abc::kWrapExisting);
        IPolyMeshSchema mesh = iPoly.getSchema();
        return (mesh.getTopologyVariance() == kHeterogenousTopology);
    }

    else if (Alembic::AbcGeom::ISubD::matches(iObj.getHeader())) {
        ISubD iSub(iObj, Alembic::Abc::kWrapExisting);
        ISubDSchema mesh = iSub.getSchema();
        return (mesh.getTopologyVariance() == kHeterogenousTopology);
    }

    return false;
}

bool isTopologyChanging(IPolyMesh poly)
{
    IPolyMeshSchema mesh = poly.getSchema();
    return (mesh.getTopologyVariance() == kHeterogenousTopology);
}

bool isTopologyChanging(std::vector<Alembic::AbcGeom::IObject> _objs)
{

	for( std::vector<Alembic::AbcGeom::IObject>::const_iterator iObj( _objs.begin() ); iObj != _objs.end(); ++iObj ) {
		if (isTopologyChanging(*iObj)) {
			return true;
		}
	}

	return false;
}

//-*****************************************************************************

bool isDeforming(IPolyMesh poly)
{
    IPolyMeshSchema mesh = poly.getSchema();
    return (mesh.getNumSamples() > 1);
}

//-*****************************************************************************

bool isNormalChanging(IPolyMesh poly)
{
    IPolyMeshSchema mesh = poly.getSchema();
    return (mesh.getNormalsParam().getNumSamples() > 1);
}

//-*****************************************************************************

bool isUVChanging(IPolyMesh poly)
{
    IPolyMeshSchema mesh = poly.getSchema();
    return (mesh.getUVsParam().getNumSamples() > 1);
}

//-*****************************************************************************

void buildBboxPrimitives(GeometryList& out, unsigned obj)
{
	// Cube vertex indices
	const unsigned VtxIndices[24] =
	{
			0, 2, 3, 1,
			0, 1, 5, 4,
			0, 4, 6, 2,
			1, 3, 7, 5,
			2, 6, 7, 3,
			4, 5, 7, 6
	};

	// Cube faces
	for (unsigned i = 0; i < 6; i++) {
		Primitive *prim = new Polygon(4, true);
		prim->vertex(3) = VtxIndices[i*4];
		prim->vertex(2) = VtxIndices[i*4+1];
		prim->vertex(1) = VtxIndices[i*4+2];
		prim->vertex(0) = VtxIndices[i*4+3];
		out.add_primitive(obj, prim);
	}
}

//-*****************************************************************************
void buildABCPrimitives(GeometryList& out, unsigned obj, const Alembic::AbcGeom::IObject iObj, chrono_t curTime)
{
	Int32ArraySamplePtr _fc;
	Int32ArraySamplePtr _fi;

	fillPrimitiveIndices(iObj, _fc, _fi, curTime);

	unsigned v_offset = 0;
	unsigned numPrimitives =_fc->size();
	// Create primitives
	for (unsigned i = 0; i < numPrimitives; i++) {
		unsigned num_verts = _fc->get()[i];
		Primitive *prim = new Polygon(num_verts, true);

		for (unsigned pv = 0; pv < num_verts; pv++) {
			prim->vertex(pv) = (*_fi)[v_offset + num_verts - pv -1]; // inverted winding order
		}

		// Add primitive to obj
		out.add_primitive(obj, prim);

		// Move offset
		v_offset += num_verts;
	}
}

void buildABCPrimitives(GeometryList& out, unsigned obj, const Alembic::AbcGeom::IPolyMesh poly, chrono_t curTime)
{
    Int32ArraySamplePtr _fc;
    Int32ArraySamplePtr _fi;

    fillPrimitiveIndices(poly, _fc, _fi, curTime);

    unsigned v_offset = 0;
    unsigned numPrimitives =_fc->size();
	// Create primitives
    for (unsigned i = 0; i < numPrimitives; i++) {
        unsigned num_verts = _fc->get()[i];
        Primitive *prim = new Polygon(num_verts, true);

        for (unsigned pv = 0; pv < num_verts; pv++) {
            prim->vertex(pv) = (*_fi)[v_offset + num_verts - pv -1]; // inverted winding order
        }

		// Add primitive to obj
        out.add_primitive(obj, prim);

		// Move offset
        v_offset += num_verts;
    }
}

//-*****************************************************************************






