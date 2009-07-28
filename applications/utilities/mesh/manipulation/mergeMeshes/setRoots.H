    argList::validArgs.clear();

    argList::validArgs.append("master root");
    argList::validArgs.append("master case");

    argList::validArgs.append("root to add");
    argList::validArgs.append("case to add");

    argList args(argc, argv);

    if (!args.check())
    {
         FatalError.exit();
    }

    fileName rootDirMaster(args.additionalArgs()[0]);
    fileName caseDirMaster(args.additionalArgs()[1]);

    fileName rootDirToAdd(args.additionalArgs()[2]);
    fileName caseDirToAdd(args.additionalArgs()[3]);

    Info<< "Master:      " << rootDirMaster << " " << caseDirMaster << nl
        << "mesh to add: " << rootDirToAdd << " " << caseDirToAdd << endl;
