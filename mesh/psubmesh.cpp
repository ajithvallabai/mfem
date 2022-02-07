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

#include <iostream>
#include <unordered_set>
#include <algorithm>
#include "psubmesh.hpp"
#include "submesh_utils.hpp"

using namespace mfem;

ParSubMesh ParSubMesh::CreateFromDomain(ParMesh &parent,
                                        Array<int> &domain_attributes)
{
   return ParSubMesh(parent, SubMesh::From::Domain, domain_attributes);
}

ParSubMesh::ParSubMesh(ParMesh &parent, SubMesh::From from,
                       Array<int> &attributes) : parent_(parent), from_(from), attributes_(attributes)
{
   MyComm = parent.GetComm();
   NRanks = parent.GetNRanks();
   MyRank = parent.GetMyRank();

   if (from == SubMesh::From::Domain)
   {
      InitMesh(parent.Dimension(), parent.SpaceDimension(), 0, 0, 0);
      auto [vtxids, elids] = SubMeshUtils::AddElementsToMesh(parent_, *this,
                                                             attributes_);
      parent_vertex_ids_ = vtxids;
      parent_element_ids_ = elids;

      // Don't let boundary conditions get generated automatically. This would
      // generate boundary elements on each rank locally, which is topologically
      // wrong for the distributed SubMesh.
      FinalizeTopology(false);
   }

   // Every rank containing elements of the SubMesh attributes now has a local
   // SubMesh. We have to connect the local meshes and assign global boundaries
   // correctly.

   DSTable v2v(parent.GetNV());
   parent.GetVertexToVertexTable(v2v);
   for (int i = 0; i < NumOfEdges; i++)
   {
      Array<int> lv;
      GetEdgeVertices(i, lv);

      // Find vertices/edge in parent mesh
      int parent_edge_id = v2v(parent_vertex_ids_[lv[0]], parent_vertex_ids_[lv[1]]);
      parent_edge_ids_.Append(parent_edge_id);
   }

   // create a GroupCommunicator on the shared edges
   GroupCommunicator sedge_comm(parent.gtopo);
   parent.GetSharedEdgeCommunicator(sedge_comm);

   Array<int> sedge_ct(sedge_comm.GroupLDofTable().Size_of_connections());

   // On each rank, locally determine if the shared edge is in the SubMesh.
   for (int i = 0; i < sedge_ct.Size(); i++)
   {
      int e1, e2;
      parent.GetFaceElements(parent.GetSharedFace(i), &e1, &e2);
      int submesh_el_id = parent_element_ids_.FindSorted(e1);
      if (submesh_el_id == -1)
      {
         // face is not in submesh
         sedge_ct[i] = 0;
      }
      else
      {
         sedge_ct[i] = 1;
      }
   }

   // Compute the sum on the root rank and broadcast the result to all ranks.
   sedge_comm.Reduce(sedge_ct, GroupCommunicator::Sum);
   sedge_comm.Bcast<int>(sedge_ct, 1);

   Array<int> pts_edge_id(parent.GetNEdges());
   pts_edge_id = -1;
   for (int i = 0; i < parent_edge_ids_.Size(); i++)
   {
      pts_edge_id[parent_edge_ids_[i]] = i;
   }

   for (int g = 1; g < parent.GetNGroups(); g++)
   {
      for (int ge = 0; ge < parent.GroupNEdges(g); ge++)
      {
         int e, o;
         parent.GroupEdge(g, ge, e, o);

         int submesh_edge = pts_edge_id[e];

         if (submesh_edge == -1)
         {
            // parent shared edge is not in SubMesh
            continue;
         }

         if (sedge_ct[ge] == 2)
         {
            // submesh_edge is shared and guaranteed interior in 2D
         }
         else if (sedge_ct[ge] == 1)
         {
            // submesh_edge is shared on one rank only in 2D, not shared anymore
         }
         else if (sedge_ct[ge] == 0)
         {
            // already handled by submesh_edge == -1
         }
      }

      // // Create global boundary
      // std::vector<int> local_bdr_faces;
      // for (int i = 0; i < faces_info.Size(); i++)
      // {
      //    if (faces_info[i].Elem2No < 0)
      //    {
      //       // this needs to map from local indices to global indices
      //       local_bdr_faces.push_back(parent_edge_ids_[i]);
      //    }
      // }

      // int recvcounts[NRanks];
      // int nlbdrfaces = (int)local_bdr_faces.size();
      // MPI_Allgather(&nlbdrfaces, 1, MPI_INT, recvcounts, 1, MPI_INT, MyComm);

      // int disp[NRanks];
      // disp[0] = 0;
      // for (int i = 1; i < NRanks; i++)
      // {
      //    disp[i] = disp[i-1] + recvcounts[i-1];
      // }

      // int ngbdrf = 0;
      // for (int i = 0; i < NRanks; i++)
      // {
      //    ngbdrf += recvcounts[i];
      // }

      // int global_bdr_faces[ngbdrf];
      // MPI_Allgatherv(local_bdr_faces.data(), (int)local_bdr_faces.size(), MPI_INT,
      //                global_bdr_faces, recvcounts, disp, MPI_INT, MyComm);

      // std::vector<std::vector<int>> lbdrfaces(NRanks);
      // for (int i = 0; i < NRanks; i++)
      // {
      //    lbdrfaces[i].resize(recvcounts[i]);
      //    for (int j = 0; j < recvcounts[i]; j++)
      //    {
      //       lbdrfaces[i][j] = global_bdr_faces[j + disp[i]];
      //    }
      //    std::sort(lbdrfaces[i].begin(), lbdrfaces[i].end());
      // }

      // for (int i = 0; i < lbdrfaces.size(); i++)
      // {
      //    if (MyRank == i) { continue; }

      //    for (int j = 0; j < lbdrfaces[MyRank].size(); j++)
      //    {
      //       if (std::binary_search(lbdrfaces[i].begin(), lbdrfaces[i].end(),
      //                              lbdrfaces[MyRank][j]))
      //       {
      //          // face was found on other processors -> not a global boundary face
      //          lbdrfaces[MyRank][j] = -1;
      //       }
      //    }
      // }

      // std::vector<int> gbdrfaces;
      // for (int i = 0; i < lbdrfaces[MyRank].size(); i++)
      // {
      //    if (lbdrfaces[MyRank][i] != -1) { gbdrfaces.push_back(lbdrfaces[MyRank][i]); }
      // }

      // out << "rank " << MyRank << ": " << gbdrfaces.size() << std::endl;

      // // The local number of (real) boundary elements on this rank
      // NumOfBdrElements = gbdrfaces.size();
      // boundary.SetSize(NumOfBdrElements);
      // be_to_edge.SetSize(NumOfBdrElements);
      // for (int i = 0, j = 0; i < gbdrfaces.size(); i++)
      // {
      //    boundary[j] = faces[gbdrfaces[i]]->Duplicate(this);
      //    be_to_edge[j++] = i;
      // }

      // FinalizeParTopo();
   }

   // Snippet that only works if shared entities are not on the boundary of the
   // submesh.
   //
   // DSTable v2v(NumOfVertices);
   // GetVertexToVertexTable(v2v);

   // Array<Connection> group_sedge_list;
   // Array<Connection> group_svert_list;
   // for (int g = 1; g < parent.GetNGroups(); g++)
   // {
   //    for (int v = 0; v < parent.GroupNVertices(g); v++)
   //    {
   //       int lvtx = parent.GroupVertex(g, v);
   //       int submesh_vtx = parent_vertex_ids_.Find(lvtx);
   //       if (submesh_vtx == -1)
   //       {
   //          continue;
   //       }
   //       svert_lvert.Append(submesh_vtx);
   //       group_svert_list.Append(Connection(g-1, svert_lvert.Size()-1));
   //    }

   //    for (int e = 0; e < parent.GroupNEdges(g); e++)
   //    {
   //       int edge, o;
   //       parent.GroupEdge(g, e, edge, o);

   //       Array<int> vert;
   //       parent.GetEdgeVertices(edge, vert);

   //       Array<int> submesh_vert(vert.Size());
   //       bool skip = false;
   //       for (int i = 0; i < vert.Size(); i++)
   //       {
   //          int submesh_vtx = parent_vertex_ids_.Find(vert[i]);
   //          if (submesh_vtx == -1)
   //          {
   //             skip = true;
   //             break;
   //          }
   //          submesh_vert[i] = submesh_vtx;
   //       }
   //       if (skip) { continue; }

   //       // mesh local edge id form vertex-to-vertex table
   //       sedge_ledge.Append(v2v(submesh_vert[0], submesh_vert[1]));
   //       shared_edges.Append(NewElement(Geometry::SEGMENT));
   //       shared_edges.Last()->SetVertices(submesh_vert);

   //       group_sedge_list.Append(Connection(g-1, shared_edges.Size()-1));
   //    }

   //    for (int e = 0; e < parent.GroupNTriangles(g); e++)
   //    {
   //       int face, o;
   //       parent.GroupTriangle(g, e, face, o);
   //       // Not relevant in 2D
   //    }

   //    for (int e = 0; e < parent.GroupNQuadrilaterals(g); e++)
   //    {
   //       int face, o;
   //       parent.GroupQuadrilateral(g, e, face, o);
   //       // Not relevant in 2D
   //    }
   // }

   // group_svert.MakeFromList(parent.GetNGroups()-1, group_svert_list);
   // group_sedge.MakeFromList(parent.GetNGroups()-1, group_sedge_list);