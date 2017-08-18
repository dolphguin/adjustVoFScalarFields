/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenFOAM Foundation
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

Application
    setWaves

Description
    Applies wave models to the entire domain for case initialisation using
    level sets for second-order accuracy.

\*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*\
Application
    adjustVoFScalarFields

Description
    Approximates a given volScalarFied using level sets for second-order
    accuracy (inspired by setWaves OpenFOAM's utility).
    It is aimed to improve the output provided by setFields when the selected
    sets do not perfectly match the goal shape e.g. when using sphereToCell in
    a fully orthogonal mesh.

Notes to myself
    Try to produce more elegant code next time.

    Finish reading/understanding the GNUGPLv3 in order to find out whether
    headers of modified OpenFOAM-files can be edited or not.

\*---------------------------------------------------------------------------*/

//#define debugON

#include "fvCFD.H"
#include "levelSet.H"
#include "pointFields.H"
#include "timeSelector.H"

#include "volPointInterpolation.H"
//#include "fixedValuePointPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions(false, false);

    Foam::argList::addOption
    (
        "alpha",
        "name",
        "name of the volume fraction field, default is \"alpha\""
    );

    #include "setRootCase.H"
    #include "createTime.H"

    instantList timeDirs = timeSelector::selectIfPresent(runTime, args);

    #include "createMesh.H"

    const pointMesh& pMesh = pointMesh::New(mesh);

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time = " << runTime.timeName() << nl << endl;

        mesh.readUpdate();

        // Read the fields which are to be set
        volScalarField alpha
        (
            IOobject
            (
                args.optionFound("alpha") ? args["alpha"] : "alpha",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ
            ),
            mesh
        );

        volScalarField alphaAfterLS
        (
            IOobject
            (
                "alphaAfterLS",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ
            ),
            mesh,
            dimensionedScalar("0", dimless, 0)
        );

        alphaAfterLS == alpha;

        pointScalarField alphaAtPoints
        (
            IOobject
            (
                "alphaAtPoints",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ
            ),
            pMesh,
            dimensionedScalar("0", dimless, 0),
            "zeroGradient"
        );

        const volPointInterpolation& vpInterp = volPointInterpolation::New(mesh);

        alphaAtPoints.ref() = vpInterp.interpolate(alpha);

        alphaAtPoints.ref() -= 0.5;

        // Set the internal fields
        alphaAfterLS.primitiveFieldRef() =
            levelSetFraction(mesh, alpha, alphaAtPoints.primitiveFieldRef(), true);

        // Set the boundary fields
        forAll(mesh.boundary(), patchi)
        {
            fvPatchScalarField& alphap = alphaAfterLS.boundaryFieldRef()[patchi];

            alphap ==
                levelSetFraction
                (
                    mesh.boundary()[patchi],
                    alpha.boundaryField()[patchi],
                    alphaAtPoints.boundaryField()[patchi].patchInternalField(),
                    true
                );
        }

        alphaAtPoints.ref() += 0.5;

        // Output
        Info<< "Writing " << alpha.name() << nl;

#ifdef debugON
        alphaAfterLS.write();
        alphaAtPoints.write();
#else
        alpha == alphaAfterLS;
        alpha.write();
#endif
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
