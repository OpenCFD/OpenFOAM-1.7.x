word fluxScheme("Kurganov");
if (mesh.schemesDict().found("fluxScheme"))
{
    fluxScheme = word(mesh.schemesDict().lookup("fluxScheme"));
    if ((fluxScheme == "Tadmor") || (fluxScheme == "Kurganov"))
    {
        Info<< "fluxScheme: " << fluxScheme << endl;
    }
    else
    {
        FatalErrorIn
        (
            "rhoCentralFoam::readFluxScheme"
        )   << "fluxScheme: " << fluxScheme
            << " is not a valid choice. "
            << "Options are: Tadmor, Kurganov"
            << abort(FatalError);
    }
}
