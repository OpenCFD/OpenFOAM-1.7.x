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
    Cloud<passiveParticle>::iterator elmnt = particles.begin();
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

forAll(lagrangianVectorNames, i)
{
    word name = lagrangianVectorNames[i];

    IOField<vector> v
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

    if (v.size())
    {
        gmvFile << name + "x" << nl;

        for (label n = 0; n < v.size(); n++)
        {
            gmvFile << v[n].x() << token::SPACE;
        }
        gmvFile << nl;

        gmvFile << name + "y" << nl;

        for (label n = 0; n < v.size(); n++)
        {
            gmvFile << v[n].y() << token::SPACE;
        }
        gmvFile << nl;

        gmvFile << name + "z" << nl;

        for (label n = 0; n < v.size(); n++)
        {
            gmvFile << v[n].z() << token::SPACE;
        }
        gmvFile << nl;
    }

}


gmvFile << "endtrace"<< nl;
