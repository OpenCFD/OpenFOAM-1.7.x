    Info<< nl << "Create Times" << endl;

    Time runTimeMaster
    (
        Time::controlDictName,
        rootDirMaster,
        caseDirMaster
    );

    Time runTimeToAdd
    (
        Time::controlDictName,
        rootDirToAdd,
        caseDirToAdd
    );
