/** @file gsMPBESBasis.cpp

    @brief instantiation file

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

#include <gsCore/gsTemplateTools.h>

#include <gsUnstructuredSplines/src/gsMPBESBasis.h>
#include <gsUnstructuredSplines/src/gsMPBESBasis.hpp>

namespace gismo
{

    // CLASS_TEMPLATE_INST gsMPBESBasis<1,real_t> ;
    CLASS_TEMPLATE_INST gsMPBESBasis<2,real_t> ;

} // end namespace gismo

