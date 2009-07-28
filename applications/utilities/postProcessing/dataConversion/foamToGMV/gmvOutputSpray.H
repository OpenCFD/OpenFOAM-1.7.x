gmvFile << "tracers " << particles.size() << nl;
for
(
    Cloud<passiveParticle>::iterator elmnt = particles.begin();
    elmnt != particles.end();
    ++elmnt
)
{
    gmvFile << elmnt().position().x() << " ";
}
gmvFile << nl;
for
(
    Cloud<passiveParticle>::iterator elmnt = particles.begin();
    elmnt != particles.end();
    ++elmnt
)
{
    gmvFile << elmnt().position().y() << " ";
}
gmvFile << nl;
for
(
    Cloud<passiveParticles>::iterator elmnt = particles.begin();
    elmnt != particles.end();
    ++elmnt
)
{
    gmvFile << elmnt().position().z() << " ";
}
gmvFile << nl;

forAll(lagrangianScalarNames, i)
{
    word name = lagrangianScalarNames[i];

    IOField<scalar> s
    (
        IOobject
        (
            name,
            runTime.timeName(),
            cloud::prefix,
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    if (s.size())
    {
        gmvFile << name  << nl;

        for (label n = 0; n < s.size(); n++)
        {
            gmvFile << s[n] << token::SPACE;
        }
        gmvFile << nl;
    }


}


gmvFile << "endtrace"<< nl;
