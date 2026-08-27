/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Description
    Template for use with codeStream.

\*---------------------------------------------------------------------------*/

#include "dictionaryEntry.H"
#include "fieldTypes.H"
#include "Ostream.H"
#include "Pstream.H"
#include "read.H"
#include "unitConversion.H"

//{{{ begin codeInclude
#line 114 "/mnt/c/Users/teppe/CFDMLProject/CFD_Project/cases/phase_1_static/system/blockMeshDict/#codeStream"

	    #include "pointField.H"
	
//}}} end codeInclude

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    void codeStream_5eb2ba567574380b33e400a92da3f469d6ab83ff
    (
        Ostream& os,
        const dictionary& dict
    )
    {
//{{{ begin code
        #line 119 "/mnt/c/Users/teppe/CFDMLProject/CFD_Project/cases/phase_1_static/system/blockMeshDict/#codeStream"

            const dictionary & d = dict;
            const scalar chord = readScalar(d.lookup("chord"));
            const scalar tOverC = readScalar(d.lookup("tOverC"));
            const scalar ffRatio = readScalar(d.lookup("farfieldRatio"));
            const scalar zSpan = readScalar(d.lookup("zSpan"));
            const scalar xMid = readScalar(d.lookup("xMid"));

	    const label n = readScalar(d.lookup("nSplinePts"));
	    


	    auto yt = [&](scalar x)
	    {
		scalar xc = x/chord;
		return 5*tOverC*chord
                (
                   0.2969*Foam::sqrt(xc)
                   - 0.1260*xc
                   - 0.3516*Foam::sqr(xc)
                   + 0.2843*Foam::pow3(xc)
                   - 0.1015*Foam::pow4(xc)
	    	);
	    };

	    // cosine-clustered stations between x0 and x1
	    auto cosineRange = [&](scalar x0, scalar x1, label n)
	    {
		scalarField x(n);

		forAll(x, i)
		{
		    scalar frac = scalar(i)/scalar(n-1);
		    scalar beta = 0.5*(1 - Foam::cos(frac*Foam::constant::mathematical::pi));
		    x[i] = x0 + beta*(x1-x0);

		}

		return x;

	    };

	    
	    auto makePts = [&](const scalarField& xs, bool upper, scalar z)
	    {
		pointField p(xs.size() - 2);
		label xsSize = xs.size();
		for(label i = 1; i < xsSize-1; ++i)
		{
		    scalar y = upper ? yt(xs[i]) : -yt(xs[i]);
		    p[i-1] = point(xs[i], y, z);
		}

		return p;

	    };

	    scalarField UF = cosineRange(0, xMid, n);
	    scalarField UA = cosineRange(xMid, chord, n);
	    scalarField LF = cosineRange(xMid, 0, n);
	    scalarField LA = cosineRange(chord, xMid, n);

	    
	    const point upperEllipse(-5*chord, 8.5745*chord, 0);
	    const point lowerEllipse(-5*chord, -8.5745*chord, 0);


	    os << "spline 0 1 " << makePts(UF, true, 0) << nl;
            os << "spline 1 2 " << makePts(UA, true, 0) << nl;
            os << "spline 2 3 " << makePts(LA, false, 0) << nl;
            os << "spline 3 0 " << makePts(LF, false, 0) << nl;

            os << "spline 12 13 " << makePts(UF, true, zSpan) << nl;
            os << "spline 13 14 " << makePts(UA, true, zSpan) << nl;
            os << "spline 14 15 " << makePts(LA, false, zSpan) << nl;
            os << "spline 15 12 " << makePts(LF, false, zSpan) << nl; 
	    
	    os << "arc 4 5 (" << upperEllipse.x() << " " << upperEllipse.y() << " 0 )" << nl;
            os << "arc 11 4 (" << lowerEllipse.x() << " " << lowerEllipse.y() << " 0 )" << nl;

            os << "arc 16 17 (" << upperEllipse.x() << " " << lowerEllipse.y() << " " << zSpan << ")" << nl;
            os << "arc 23 16 (" << lowerEllipse.x() << " " << lowerEllipse.y() << " " << zSpan << ")" << nl;
	

	
//}}} end code
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

