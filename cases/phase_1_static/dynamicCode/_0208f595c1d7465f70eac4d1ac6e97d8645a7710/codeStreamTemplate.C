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
#line 38 "/mnt/c/Users/teppe/CFDMLProject/CFD_Project/cases/phase_1_static/system/blockMeshDict/#codeStream"

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
    void codeStream_0208f595c1d7465f70eac4d1ac6e97d8645a7710
    (
        Ostream& os,
        const dictionary& dict
    )
    {
//{{{ begin code
        #line 43 "/mnt/c/Users/teppe/CFDMLProject/CFD_Project/cases/phase_1_static/system/blockMeshDict/#codeStream"

	    const dictionary & d = dict;
  	    const scalar chord = readScalar(d.lookup("chord"));
	    const scalar tOverC = readScalar(d.lookup("tOverC"));
	    const scalar ffRatio = readScalar(d.lookup("farfieldRatio"));
	    const scalar zSpan = readScalar(d.lookup("zSpan"));
	    const scalar xMid = readScalar(d.lookup("xMid"));	    

	
	    auto yt = [&](scalar x)
	    {
		scalar xc = x/chord;
	 	return 5*tOverC*chord*
		(
		   0.2969*Foam::sqrt(xc)
		   - 0.1260*xc
		   - 0.3516*Foam::sqr(xc)
		   + 0.2843*Foam::pow3(xc)
		   - 0.1015*Foam::pow4(xc)

		);

	    };


	   point WE(0, 0, 0);
	   point NE(xMid, yt(xMid), 0);
	   point EE(chord, 0, 0);
  	   point SE(xMid, -yt(xMid), 0);

	   point oWE(-10*chord, 0, 0);
           point oNE1(xMid, 10*chord, 0);
           point oNE2(chord, 10*chord, 0);
           point oNE3(21*chord, 10*chord, 0);
	   point oSE1(xMid, -10*chord, 0);
	   point oSE2(chord, -10*chord, 0);
           point oSE3(21*chord, -10*chord, 0);
           point oEE(21*chord, 0, 0);

	   pointField pts(24);
	   pointField layer0({WE, NE, EE, SE, oWE, oNE1, oNE2, oNE3,
			      oEE, oSE3, oSE2, oSE1});

	   forAll(layer0, i)
	   {
	      pts[i] = layer0[i];
	      pts[i+12] = layer0[i] + point(0, 0, zSpan);
	   }

	   forAll(pts, i)
	   {
		os << pts[i] << nl;

	   }

	
//}}} end code
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

