/*
 * TinyEKF: Extended Kalman Filter for Arduino and TeensyBoard.
 *
 * Copyright (C) 2015 Simon D. Levy
 *
 * This code is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This code is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this code.  If not, see <http:#www.gnu.org/licenses/>.
 */

#include "TinyEKF.h"

TinyEKF::TinyEKF(int n, int m) 
{
    ekf_init(&this->ekf, n, m);

    this->x = this->ekf.x;
}

TinyEKF::~TinyEKF()
{
}

void TinyEKF::setP(int i, int j, double value)
{
    ekf_setP(&this->ekf, i, j, value);
}

void TinyEKF::setQ(int i, int j, double value)
{
    ekf_setQ(&this->ekf, i, j, value);
}

void TinyEKF::setR(int i, int j, double value)
{
    ekf_setR(&this->ekf, i, j, value);
}

double TinyEKF::getX(int i)
{
    return this->ekf.x[i];
}

bool TinyEKF::step(double * z)
{        
    this->model(this->ekf.fx, this->ekf.F, this->ekf.hx, this->ekf.H);

    return ekf_step(&this->ekf, z) == 0 ? true : false;
}

void TinyEKF::set(double * A, int i, int j, double value)
{
    ekf_set(&this->ekf, A, i, j, value);
}

