/*
Copyright (C) 2014, 2015, 2016
Rafael Guglielmetti, rafael.guglielmetti@unifr.ch
*/

/*
This file is part of AlVin.

CoxIter is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

CoxIter is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with AlVin. If not, see <http://www.gnu.org/licenses/>.
*/

/*!
 * \file rationalinteger_notreflective.h
 * \author Rafael Guglielmetti
 *
 * \class RationalInteger_NotReflective
 * \brief Try to create systems of equations to show that the rational quadratic
 * form is not reflective
 */

#ifndef RATIONALINTEGER_NOTREFLECTIVE_H
#define RATIONALINTEGER_NOTREFLECTIVE_H

#include "notreflective.h"

class RationalInteger_NotReflective : public NotReflective {
public:
  RationalInteger_NotReflective(AlVin *v);

  void createSystemEquations(NotReflective_Graph nrg);
};

#endif // RATIONALINTEGER_NOTREFLECTIVE_H
