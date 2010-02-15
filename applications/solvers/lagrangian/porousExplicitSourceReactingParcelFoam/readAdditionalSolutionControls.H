dictionary additional = mesh.solutionDict().subDict("additional");

bool dpdt = additional.lookupOrDefault("dpdt", true);
bool eWork = additional.lookupOrDefault("eWork", true);
bool hWork = additional.lookupOrDefault("hWork", true);
