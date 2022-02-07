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

#include "submesh_utils.hpp"

namespace mfem
{
namespace SubMeshUtils
{
int UniqueIndexGenerator::Get(int i, bool &new_index)
{
   auto f = idx.find(i);
   if (f == idx.end())
   {
      idx[i] = counter;
      new_index = true;
      return counter++;
   }
   else
   {
      new_index = false;
      return (*f).second;
   }
}

bool IsSubMesh(const Mesh *m)
{
   return dynamic_cast<const SubMesh *>(m) != nullptr;
}

bool ElementHasAttribute(const Element &el, const Array<int> &attributes)
{
   for (int a = 0; a < attributes.Size(); a++)
   {
      if (el.GetAttribute() == attributes[a])
      {
         return true;
      }
   }
   return false;
}

std::tuple<Array<int>, Array<int>>
                                AddElementsToMesh(const Mesh& parent, Mesh& mesh, const Array<int> &attributes)
{
   Array<int> parent_vertex_ids, parent_element_ids;
   SubMeshUtils::UniqueIndexGenerator vertex_ids;
   for (int i = 0; i < parent.GetNE(); i++)
   {
      const Element *pel = parent.GetElement(i);

      if (!SubMeshUtils::ElementHasAttribute(*pel, attributes)) { continue; }

      Array<int> v;
      pel->GetVertices(v);
      Array<int> submesh_v(v.Size());

      for (int iv = 0; iv < v.Size(); iv++)
      {
         bool new_vertex;
         int mesh_vertex_id = v[iv];
         int submesh_vertex_id = vertex_ids.Get(mesh_vertex_id, new_vertex);
         if (new_vertex)
         {
            mesh.AddVertex(parent.GetVertex(mesh_vertex_id));
            parent_vertex_ids.Append(mesh_vertex_id);
         }
         submesh_v[iv] = submesh_vertex_id;
      }

      Element *el = mesh.NewElement(parent.GetElementType(i));
      el->SetVertices(submesh_v);
      el->SetAttribute(pel->GetAttribute());
      mesh.AddElement(el);
      parent_element_ids.Append(i);
   }
   return std::tuple<Array<int>, Array<int>>(parent_vertex_ids,
                                             parent_element_ids);
}

};
};
