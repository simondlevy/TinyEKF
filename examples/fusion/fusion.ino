#include "TinyEKF.h"

class Fuser : public TinyEKF {

    public:

        // Eight state values, four measurement values
        Fuser() : TinyEKF(1, 1)
    {            
        this->setP(0, 0, .01);
        this->setQ(0, 0, .01);
        this->setR(0, 0, .01);

    }

    protected:

        void model(float * fx, float * F, float * hx, float * H)
        {
            fx[0] = this->x[0];
            hx[0] = fx[0];

            set(F, 0, 0, 1);
            set(H, 0, 0, 1);
        }
};

Fuser ekf;

void setup() {

    Serial.begin(9600);
}


void loop() {

    static int count;
    const int LOOPSIZE = 1000;

    float z[1];

    z[0] = sin(2*M_PI*count/LOOPSIZE);

    ekf.step(z);

    Serial.println(ekf.getX(0));

    count = (count + 1) % LOOPSIZE;
}
