// Copyright (c) 2010-2021, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.


#ifndef MFEM_SUBMESH_UTILS
#define MFEM_SUBMESH_UTILS

#include <unordered_map>
#include "submesh.hpp"
#include "psubmesh.hpp"

namespace mfem
{
// TODO: Decide if these are just "submesh" utilities or really general "mesh"
// utilities
namespace SubMeshUtils
{
struct UniqueIndexGenerator
{
    int counter = 0;
    std::unordered_map<int, int> idx;
    int Get(int i, bool &new_index);
};

bool IsSubMesh(const Mesh *m);

bool ElementHasAttribute(const Element &el, const Array<int> &attributes);

std::tuple<Array<int>, Array<int>> AddElementsToMesh(const Mesh& parent, Mesh& mesh, const Array<int> &attributes);
};
};

#endif
