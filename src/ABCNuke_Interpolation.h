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

#ifndef _ABCNuke_Interpolation_h_
#define _ABCNuke_Interpolation_h_

#include "DDImage/Vector3.h"
#include "DDImage/Matrix4.h"
#include "DDImage/Quaternion.h"
#include "DDImage/DDMath.h"
#include "ABCNuke_MatrixHelper.h"

using namespace DD::Image;

double getWeightAndIndex(double iFrame,
    Alembic::AbcCoreAbstract::TimeSamplingPtr iTime, size_t numSamps,
    Alembic::AbcCoreAbstract::index_t & oIndex,
    Alembic::AbcCoreAbstract::index_t & oCeilIndex);

// Not needed. Using OpenEXR::Imath functions for now
DD::Image::Matrix4 Matrix4_lerp(const DD::Image::Matrix4 &start_mtx, const DD::Image::Matrix4 &end_mtx, double amt);

template<typename T>
        T simpleLerp(double alpha, T val1, T val2)
{
    double dv = static_cast<double>( val1 );
    return static_cast<T>( dv + alpha * (static_cast<double>(val2) - dv) );
}
#endif
