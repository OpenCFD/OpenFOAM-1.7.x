    argList::validArgs.clear();

    argList::validArgs.append("master root");
    argList::validArgs.append("master case");
    argList::validOptions.insert("masterRegion", "name");

    argList::validArgs.append("root to add");
    argList::validArgs.append("case to add");
    argList::validOptions.insert("addRegion", "name");

    argList args(argc, argv);

    if (!args.check())
    {
         FatalError.exit();
    }

    fileName rootDirMaster(args.additionalArgs()[0]);
    fileName caseDirMaster(args.additionalArgs()[1]);
    word masterRegion = polyMesh::defaultRegion;
    args.optionReadIfPresent("masterRegion", masterRegion);

    fileName rootDirToAdd(args.additionalArgs()[2]);
    fileName caseDirToAdd(args.additionalArgs()[3]);
    word addRegion = polyMesh::defaultRegion;
    args.optionReadIfPresent("addRegion", addRegion);

    Info<< "Master:      " << rootDirMaster << " " << caseDirMaster
        << "  region " << masterRegion << nl
        << "mesh to add: " << rootDirToAdd << " " << caseDirToAdd
        << "  region " << addRegion << endl;
