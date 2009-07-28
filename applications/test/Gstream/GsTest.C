#include "Xstream.H"
#include "GLstream.H"
#include "PSstream.H"
#include "shapes2D.H"
#include "IStringStream.H"

using namespace Foam;

int main()
{
    colour mauve("mauve", 1, 0, 1);

    lineStyle solid("Solid", 2.0, 2.0, IStringStream("1(1.0)")());
    lineStyle broken("Broken", 2.0, 10.0, IStringStream("4(1 1 4 1)")());


    //Xstream wind
    GLstream wind
    //PSstream wind
    (
        "GsTest",
        primary("Black"),
        primary("White"),
        0.5, 0.5, 0.5, 0.5, 500, 500
    );

    do
    {
        wind << rectangle2D(point2D(0.0, 0.0), point2D(100.0, 100.0));

        wind.setColour(mauve);
        wind.setLineStyle(solid);

        wind << line2D(point2D(0.0, 0.0), point2D(0.0, 200.0));
        wind << line2D(point2D(0.0, 200.0), point2D(200.0, 200.0));

        //wind.setLineStyle(broken);

        wind << line2D(point2D(200.0, 200.0), point2D(200.0, 0.0));
        wind << line2D(point2D(200.0, 0.0), point2D(0.0, 0.0));

        wind << string2D(point2D(200.0, 0.0), "Hi there");

    } while (wind.waitForEvent());

    return 0;
}
