{

bool foundCell = false;
if (injectorCell >= 0)
{
    const vector& C = mesh_.C()[injectorCell];
    injectionPosition += 1.0e-9*(C - injectionPosition);
    
    foundCell = mesh_.pointInCell
    (
        injectionPosition,
        injectorCell
    );
}

reduce(foundCell, orOp<bool>());

// NN.
// Due to the disc injection, the injectionPosition
// can be moved outside the domain.
// try point injection

if (!foundCell)
{
    injectionPosition = it->position(n);
    injectorCell = mesh_.findCell(injectionPosition);
    
    if (injectorCell >= 0)
    {
        const vector& C = mesh_.C()[injectorCell];
        injectionPosition += 1.0e-6*(C - injectionPosition);

        foundCell = mesh_.pointInCell
        (
            injectionPosition,
            injectorCell
        );
    }
    reduce(foundCell, orOp<bool>());

    // if point injection also failed then
    // find nearest cell and try that one
    // the point is probably on an edge
    if (!foundCell)
    {
        injectorCell =
            mesh_.findNearestCell(injectionPosition);
        
        if (injectorCell >= 0)
        {
            
            const vector& C = mesh_.C()[injectorCell];
            injectionPosition += 1.0e-9*(C - injectionPosition);
            
            foundCell = mesh_.pointInCell
            (
                injectionPosition,
                injectorCell
            );
        }
        reduce(foundCell, orOp<bool>());
        
        if (!foundCell)
        {
            FatalError 
                << "Cannot find injection position " 
                    << injectionPosition
                    << abort(FatalError);
        }
    }
}

}
