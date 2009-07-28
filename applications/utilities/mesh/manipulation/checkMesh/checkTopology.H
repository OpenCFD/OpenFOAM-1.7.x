#include "label.H"
#include "wordList.H"

namespace Foam
{
    class polyMesh;

    bool checkSync(const wordList& names);

    label checkTopology(const polyMesh&, const bool, const bool);
}
