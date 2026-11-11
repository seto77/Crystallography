namespace Crystallography;

//汎用性の高い列挙体をここで定義

public enum LengthUnitEnum
{
    None,

    Meter,
    CentiMeter,
    MilliMeter,
    MicroMeter,
    NanoMeter,
    Angstrom,
    PicoMeter,

    MeterInverse,
    CentiMeterInverse,
    MilliMeterInverse,
    MicroMeterInverse,
    NanoMeterInverse,
    AngstromInverse,
    PicoMeterInverse
}

public enum AngleUnitEnum
{
    Degree,
    Radian,
    CentiDegree,
    MilliRadian
}

public enum FourierDirectionEnum
{
    Forward,
    Inverse
}

public enum VH_DirectionEnum
{
    Vertical,
    Horizontal
}

public enum EnergyUnitEnum
{
    eV,
    KeV,
    MeV
}

public enum TimeUnitEnum
{
    Seccond,
    MilliSecond,
    MicroSecond,
    NanoSecond
}

// 261111Cl 変更前: 中身の無いクラスが名前空間 Crystallography の中で System.Enum を隠していた (参照は 0 件。System.Enum は完全修飾で使われている)
//internal class Enum
//{
//}