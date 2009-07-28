HashSet<word> volScalarHash;
HashSet<word> volVectorHash;
HashSet<word> surfScalarHash;
HashSet<word> surfVectorHash;
HashSet<word> sprayScalarHash;
HashSet<word> sprayVectorHash;

forAll(timeDirs, timeI)
{
    runTime.setTime(timeDirs[timeI], timeI);

    // Add all fields to hashtable
    IOobjectList objects(mesh, runTime.timeName());
    {
        wordList fieldNames(objects.names(volScalarField::typeName));
        forAll(fieldNames, fieldI)
        {
            volScalarHash.insert(fieldNames[fieldI]);
        }
    }
    {
        wordList fieldNames(objects.names(volVectorField::typeName));
        forAll(fieldNames, fieldI)
        {
            volVectorHash.insert(fieldNames[fieldI]);
        }
    }
    {
        wordList fieldNames(objects.names(surfaceScalarField::typeName));
        forAll(fieldNames, fieldI)
        {
            surfScalarHash.insert(fieldNames[fieldI]);
        }
    }
    {
        wordList fieldNames(objects.names(surfaceVectorField::typeName));
        forAll(fieldNames, fieldI)
        {
            surfVectorHash.insert(fieldNames[fieldI]);
        }
    }


    // Same for spray
    IOobjectList sprayObjects(mesh, runTime.timeName(), cloud::prefix);
    {
        wordList fieldNames(sprayObjects.names(scalarIOField::typeName));
        forAll(fieldNames, fieldI)
        {
            sprayScalarHash.insert(fieldNames[fieldI]);
        }
    }
    {
        wordList fieldNames(sprayObjects.names(vectorIOField::typeName));
        forAll(fieldNames, fieldI)
        {
            sprayVectorHash.insert(fieldNames[fieldI]);
        }
    }
}


wordList volScalarNames(volScalarHash.toc());
wordList volVectorNames(volVectorHash.toc());
wordList surfScalarNames(surfScalarHash.toc());
wordList surfVectorNames(surfVectorHash.toc());
wordList sprayScalarNames(sprayScalarHash.toc());
wordList sprayVectorNames(sprayVectorHash.toc());
