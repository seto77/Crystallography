using System;
using System.Drawing;

namespace Crystallography
{
    public class PhotoInformation
    {
        public double L1 = 0, L2 = 0, L3 = 0, Theta = 0;
        public double Tilt1 = 0, Tilt2 = 0, Tilt1Err = 0, Tilt2Err = 0;

        public double d1 = 0, d2 = 0, d3 = 0;
        public double d1max, d1min, d2max, d2min, d3max, d3min, theta_err;
        public bool IsTriangleMode;
        public PointF P1, P2;
        public bool Paintable;

        public PhotoInformation(double l1, double l2, double l3, double theta, double l1err, double l2err, double l3err, double theta_err, double tilt1, double tilt2, double tilt1err, double tilt2err,
            bool isTriangleMode, double waveLength, double cameraLength)
        {
            Paintable = true;
            L1 = l1;
            L2 = l2;
            L3 = l3;
            Theta = theta;
            this.theta_err = theta_err;
            Tilt1 = tilt1;
            Tilt2 = tilt2;
            Tilt1Err = tilt1err;
            Tilt2Err = tilt2err;

            IsTriangleMode = isTriangleMode;
            if (IsTriangleMode)//三辺モードのとき
            {
                if (L3 < L1 + L2 && L2 < L1 + L3 && L1 < L2 + L3)
                    Theta = Math.Acos((L1 * L1 + L2 * L2 - L3 * L3) / 2 / L1 / L2);
                else
                {
                    Theta = 0;
                    Paintable = false;
                }
            }
            else//二辺挟角モードのとき
            {
                if (0 < theta && theta < Math.PI)
                    L3 = Math.Sqrt(L1 * L1 + L2 * L2 - 2 * L1 * L2 * Math.Cos(Theta));
                else
                {
                    L3 = 0;
                    Paintable = false;
                }
            }

            if (Paintable)//描画可能であれば点の位置を決定
            {
                P1 = new PointF((float)L1, 0f);
                P2 = new PointF((float)(L2 * Math.Cos(Theta)), (float)(-L2 * Math.Sin(Theta)));

                d1 = waveLength / 2 / Math.Sin(Math.Atan(L1 / cameraLength) / 2);
                d1max = waveLength / 2 / Math.Sin(Math.Atan(L1 * (1 - l1err) / cameraLength) / 2);
                d1min = waveLength / 2 / Math.Sin(Math.Atan(L1 * (1 + l1err) / cameraLength) / 2);

                d2 = waveLength / 2 / Math.Sin(Math.Atan(L2 / cameraLength) / 2);
                d2max = waveLength / 2 / Math.Sin(Math.Atan(L2 * (1 - l2err) / cameraLength) / 2);
                d2min = waveLength / 2 / Math.Sin(Math.Atan(L2 * (1 + l2err) / cameraLength) / 2);

                d3 = waveLength / 2 / Math.Sin(Math.Atan(L3 / cameraLength) / 2);
                d3max = waveLength / 2 / Math.Sin(Math.Atan(L3 * (1 - l3err) / cameraLength) / 2);
                d3min = waveLength / 2 / Math.Sin(Math.Atan(L3 * (1 + l3err) / cameraLength) / 2);
            }
        }

        public void Rot(double angle)
        {
            P1 = new PointF((float)(Math.Cos(angle) * P1.X + Math.Sin(angle) * P1.Y), (float)(Math.Sin(angle) * P1.X + Math.Cos(angle) * P1.Y));
            P2 = new PointF((float)(Math.Cos(angle) * P2.X + Math.Sin(angle) * P2.Y), (float)(Math.Sin(angle) * P2.X + Math.Cos(angle) * P2.Y));
        }
    }
}