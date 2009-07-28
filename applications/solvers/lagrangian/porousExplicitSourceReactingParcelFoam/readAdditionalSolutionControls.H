dictionary additional = mesh.solutionDict().subDict("additional");

bool dpdt = true;
if (additional.found("dpdt"))
{
    additional.lookup("dpdt") >> dpdt;
}

bool eWork = true;
if (additional.found("eWork"))
{
    additional.lookup("eWork") >> eWork;
}

bool hWork = true;
if (additional.found("hWork"))
{
    additional.lookup("hWork") >> hWork;
}

