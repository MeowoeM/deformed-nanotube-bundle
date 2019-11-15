module ClassBond
use class_kBase

    type Bond
        !!<member name="indexTube1">The index of one nanotube.</member>
        integer, private :: indexTube1
        !!<member name="indexTube2">The index of one nanotube.</member>
        integer, private :: indexTube2
        !!<member name="base1">k base for tube 1</member>
        type(kBase), private :: base1
        !!<member name="base2">k base for tube 2</member>
        type(kBase), private :: base2
        !!<member name="angleFromX">the angle from x direction to the line connecting the centers of two tubes.</member>
        doubleprecision :: angleFromX
    end type

    interface Bond
        procedure :: constructBond
    end interface
end module
