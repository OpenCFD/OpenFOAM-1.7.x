#include "label.H"
#include "HashSet.H"
#include "labelVector.H"

namespace Foam
{
    class polyMesh;
    class wedgePolyPatch;

    label findOppositeWedge(const polyMesh&, const wedgePolyPatch&);

    bool checkWedges
    (
        const polyMesh&,
        const bool report,
        const Vector<label>&,
        labelHashSet*
    );

    label checkGeometry(const polyMesh& mesh, const bool allGeometry);
}
