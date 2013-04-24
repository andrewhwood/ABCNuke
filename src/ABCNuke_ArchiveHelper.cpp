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
#include "ABCNuke_ArchiveHelper.h"
//-*****************************************************************************
#include <stdio.h>
#include <iostream>

using namespace Alembic::AbcGeom;

//kh
void getABCGeosXforms(Alembic::Abc::IObject & iObj,
                      std::vector<Alembic::AbcGeom::IObject> & _objs,
                      std::vector<Alembic::AbcGeom::IPolyMesh> & _meshes,
                      std::vector< std::vector<Alembic::AbcGeom::IXform> > & _xforms)
{


    unsigned int numChildren = iObj.getNumChildren();

    for (unsigned i=0; i<numChildren; ++i)
    {
        IObject child( iObj.getChild( i ));
        if ( Alembic::AbcGeom::IPolyMesh::matches(child.getHeader()))
        //     || Alembic::AbcGeom::ISubD::matches(child.getHeader())) 
        {
            // add object
            _objs.push_back(child);
            
            // add mesh
            IPolyMesh poly(child, Alembic::Abc::kWrapExisting);
            _meshes.push_back(poly);

            // add all the xforms for the object to the list
            std::vector<Alembic::AbcGeom::IXform> xfs;
            
            IObject parent = iObj;

            while ( parent )
            {
                if ( IXform::matches( parent.getHeader() ) )
                {
                    IXform xf( parent, kWrapExisting );
                    xfs.push_back( xf );
                }
                parent = parent.getParent();

            }
            _xforms.push_back( xfs );
        }

        
        if (child.getNumChildren() > 0) 
        {
            getABCGeosXforms(child, _objs, _meshes, _xforms);
        }
        
    }
}
// end kh

void getABCGeos(Alembic::Abc::IObject & iObj,
		std::vector<Alembic::AbcGeom::IObject> & _objs)
{

	unsigned int numChildren = iObj.getNumChildren();

	for (unsigned i=0; i<numChildren; ++i)
	{
		IObject child( iObj.getChild( i ));
		if ( Alembic::AbcGeom::IPolyMesh::matches(child.getHeader())
		|| Alembic::AbcGeom::ISubD::matches(child.getHeader())) {
			_objs.push_back(child);
		}

		if (child.getNumChildren() > 0) {
			getABCGeos(child, _objs);
		}
	}
}

//-*****************************************************************************

void getABCXforms(Alembic::Abc::IObject & iObj,
		std::vector<Alembic::AbcGeom::IXform> & _objs)
{

	unsigned int numChildren = iObj.getNumChildren();

	for (unsigned i=0; i<numChildren; ++i)
	{
		IObject child( iObj.getChild( i ));
		if ( Alembic::AbcGeom::IXform::matches(child.getHeader()) ) {
			IXform xf(child, Alembic::Abc::kWrapExisting);
			_objs.push_back(xf);
		}

		if (child.getNumChildren() > 0) {
			getABCXforms(child, _objs);
		}
	}
}

//-*****************************************************************************

void getABCCameras(Alembic::Abc::IObject & iObj,
		std::vector<Alembic::AbcGeom::ICamera> & _objs)
{

	unsigned int numChildren = iObj.getNumChildren();

	for (unsigned i=0; i<numChildren; ++i)
	{
		IObject child( iObj.getChild( i ));
		if ( Alembic::AbcGeom::ICamera::matches(child.getHeader()) ) {
			ICamera cam(child, Alembic::Abc::kWrapExisting);
			_objs.push_back(cam);
		}

		if (child.getNumChildren() > 0) {
			getABCCameras(child, _objs);
		}
	}
}
//-*****************************************************************************

bool getNamedXform( IObject iObjTop, const std::string &iName, IXform &iXf )
{
	// Return true if found

	const Alembic::AbcGeom::MetaData &md = iObjTop.getMetaData();
	if ( (iObjTop.getName() == iName) && (IXform::matches( md )) )
	{
		iXf = IXform(iObjTop, kWrapExisting );
		return true;
	}

	// now the child objects
	for ( size_t i = 0 ; i < iObjTop.getNumChildren() ; i++ )
	{
		if (getNamedXform(IObject( iObjTop, iObjTop.getChildHeader( i ).getName() ), iName, iXf ))
			return true;
	}

	return false;

}

//-*****************************************************************************

bool getNamedCamera( IObject iObjTop, const std::string &iName, ICamera &iCam )
{
	// Return true if found

	const Alembic::AbcGeom::MetaData &md = iObjTop.getMetaData();
	if ( (iObjTop.getName() == iName) && (ICamera::matches( md )) )
	{
		iCam = ICamera(iObjTop, kWrapExisting );
		return true;
	}

	// now the child objects
	for ( size_t i = 0 ; i < iObjTop.getNumChildren() ; i++ )
	{
		if (getNamedCamera(IObject( iObjTop, iObjTop.getChildHeader( i ).getName() ), iName, iCam ))
			return true;
	}

	return false;

}
//-*****************************************************************************

void getXformTimeSpan(IXform iXf, chrono_t& first, chrono_t& last, bool inherits) {

	IXformSchema xf = iXf.getSchema();
	TimeSamplingPtr ts = xf.getTimeSampling();
	first = std::min(first, ts->getSampleTime(0) );
	last = std::max(last, ts->getSampleTime(xf.getNumSamples()-1) );
	if (inherits && xf.getInheritsXforms()) {
		IObject parent = iXf.getParent();

		// Once the Archive's Top Object is reached, IObject::getParent() will
		// return an invalid IObject, and that will evaluate to False.
		while ( parent )
		{
			if ( Alembic::AbcGeom::IXform::matches(parent.getHeader()) ) {
				IXform x( parent, kWrapExisting );
				getXformTimeSpan(x, first, last, inherits);
			}
		}

	}
}

//-*****************************************************************************

void getCameraTimeSpan(ICamera iCam, chrono_t& first, chrono_t& last) {

	ICameraSchema cam = iCam.getSchema();
	TimeSamplingPtr ts = cam.getTimeSampling();
	first = std::min(first, ts->getSampleTime(0) );
	last = std::max(last, ts->getSampleTime(cam.getNumSamples()-1) );
}

//-*****************************************************************************

void getPolyMeshTimeSpan(IPolyMesh iPoly, chrono_t& first, chrono_t& last) {

	IPolyMeshSchema mesh = iPoly.getSchema();
	TimeSamplingPtr ts = mesh.getTimeSampling();
	first = std::min(first, ts->getSampleTime(0) );
	last = std::max(last, ts->getSampleTime(mesh.getNumSamples()-1) );
}

//-*****************************************************************************

void getSubDTimeSpan(ISubD iSub, chrono_t& first, chrono_t& last) {

	ISubDSchema mesh = iSub.getSchema();
	TimeSamplingPtr ts = mesh.getTimeSampling();
	first = std::min(first, ts->getSampleTime(0) );
	last = std::max(last, ts->getSampleTime(mesh.getNumSamples()-1) );
}

//-*****************************************************************************

void getObjectTimeSpan(IObject obj, chrono_t& first, chrono_t& last, bool doChildren)
{
	if ( Alembic::AbcGeom::IPolyMesh::matches(obj.getHeader()) ) {
		IPolyMesh iPoly(obj, Alembic::Abc::kWrapExisting);
		getPolyMeshTimeSpan(iPoly, first, last);
	}

	else if ( Alembic::AbcGeom::ISubD::matches(obj.getHeader()) ) {
		ISubD iSub(obj, Alembic::Abc::kWrapExisting);
		getSubDTimeSpan(iSub, first, last);
	}

	else if ( Alembic::AbcGeom::IXform::matches(obj.getHeader()) ) {
		IXform iXf(obj, Alembic::Abc::kWrapExisting);
		getXformTimeSpan(iXf, first, last, false);
	}

	else if ( Alembic::AbcGeom::ICamera::matches(obj.getHeader()) ) {
		ICamera iCam(obj, Alembic::Abc::kWrapExisting);
		getCameraTimeSpan(iCam, first, last);
	}

	if (doChildren) {
		// do this object's children too
		for (unsigned i=0; i < obj.getNumChildren(); ++i)
		{
			IObject child( obj.getChild( i ));
			getObjectTimeSpan(child, first, last, doChildren);
		}
	}
}
//-*****************************************************************************
void getABCTimeSpan(IArchive archive, chrono_t& first, chrono_t& last)
{
	// TO DO: Is the childBounds property reliable to get the full archive's span?

	if (!archive.valid())
		return;

	IObject archiveTop = archive.getTop();

	unsigned int numChildren = archiveTop.getNumChildren();

	for (unsigned i=0; i<numChildren; ++i)
	{
		IObject obj( archiveTop.getChild( i ));
		getObjectTimeSpan(obj, first, last, true);

	}

}
//-*****************************************************************************
