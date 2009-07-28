{
    /* Output the magic number. */
    writeInt(fvParticleFile, FV_MAGIC);

    /* Output file header and version number. */
    writeStr80(fvParticleFile, "FVPARTICLES");

    /* version */
    writeInt(fvParticleFile, 1);
    writeInt(fvParticleFile, 1);

    int nFields = sprayScalarNames.size() + 3*sprayVectorNames.size();
    writeInt(fvParticleFile, nFields);

    forAll(sprayScalarNames, nameI)
    {
        writeStr80(fvParticleFile, sprayScalarNames[nameI]);
    }
    forAll(sprayVectorNames, nameI)
    {
        const string& name(sprayVectorNames[nameI]);
        writeStr80(fvParticleFile, name + ("x;" + name));
        writeStr80(fvParticleFile, name + ("y;" + name));
        writeStr80(fvParticleFile, name + ("z;" + name));
    }
}
