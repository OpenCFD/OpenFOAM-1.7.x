gmvFile << "tracers " << particles.size() << nl;
for
(
    discretePhase::iterator elmnt = particles.begin();
    elmnt != particles.end();
    ++elmnt
)
{
    gmvFile << elmnt().position().x() << " ";
}
gmvFile << nl;
for
(
    discretePhase::iterator elmnt = particles.begin();
    elmnt != particles.end();
    ++elmnt
)
{
    gmvFile << elmnt().position().y() << " ";
}
gmvFile << nl;
for
(
    discretePhase::iterator elmnt = particles.begin();
    elmnt != particles.end();
    ++elmnt
)
{
    gmvFile << elmnt().position().z() << " ";
}
gmvFile << nl;

gmvFile << "U" << nl;
for
(
    discretePhase::iterator elmnt = particles.begin();
    elmnt != particles.end();
    ++elmnt
)
{
    gmvFile << elmnt().velocity().x() << " ";
}
gmvFile << nl;
gmvFile << "V" << nl;
for
(
    discretePhase::iterator elmnt = particles.begin();
    elmnt != particles.end();
    ++elmnt
)
{
    gmvFile << elmnt().velocity().y() << " ";
}
gmvFile << nl;

gmvFile << "W" << nl;
for
(
    discretePhase::iterator elmnt = particles.begin();
    elmnt != particles.end();
    ++elmnt
)
{
    gmvFile << elmnt().velocity().z() << " ";
}
gmvFile << nl;

gmvFile << "Diam" << nl;

for
(
    discretePhase::iterator elmnt = particles.begin();
    elmnt != particles.end();
    ++elmnt
)
{
    gmvFile << elmnt().d() << " ";
}


gmvFile << "endtrace"<< nl;
