    IOdictionary pdfDictionary
    (
        IOobject
        (
            "pdfDictionary",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    label nIntervals
    (
        readLabel(pdfDictionary.lookup("nIntervals"))
    );

    label nSamples
    (
        readLabel(pdfDictionary.lookup("nSamples"))
    );

    label samples[nIntervals];

    for(label i=0;i<nIntervals;i++)
    {
        samples[i] = 0;
    }

