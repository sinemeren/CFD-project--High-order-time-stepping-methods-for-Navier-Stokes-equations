#include "helper.hpp"
#include "visual.hpp"
#include <stdio.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkStructuredGrid.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkTuple.h>
#include <vtkStructuredGridWriter.h>
#include "vector"
#include <bits/stdc++.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <iostream>

void write_vtkFile(const char *szProblem,
		 int    timeStepNumber,
		 double xlength,
                 double ylength,
                 int    imax,
                 int    jmax,
		 double dx,
		 double dy,
                 matrix<double> &U,
                 matrix<double> &V,
                 matrix<double> &P) {

        int i,j;
  char szFileName[80];
  FILE *fp=NULL;
  sprintf( szFileName, "%s.%i.vtk", szProblem, timeStepNumber );
  fp = fopen( szFileName, "w");
  if( fp == NULL )		       
  {
    char szBuff[80];
    sprintf( szBuff, "Failed to open %s", szFileName );
    ERROR( szBuff );
    return;
  }

  write_vtkHeader( fp, imax, jmax, dx, dy);
  write_vtkPointCoordinates(fp, imax, jmax, dx, dy);

  fprintf(fp,"POINT_DATA %i \n", (imax+1)*(jmax+1) );
	
  fprintf(fp,"\n");
  fprintf(fp, "VECTORS velocity float\n");
  for(j = 0; j < jmax+1; j++) {
    for(i = 0; i < imax+1; i++) {
      fprintf(fp, "%f %f 0\n", (U[i][j] + U[i][j+1]) * 0.5, (V[i][j] + V[i+1][j]) * 0.5 );
    }
  }

  fprintf(fp,"\n");
  fprintf(fp,"CELL_DATA %i \n", ((imax)*(jmax)) );
  fprintf(fp, "SCALARS pressure float 1 \n");
  fprintf(fp, "LOOKUP_TABLE default \n");
  for(j = 1; j < jmax+1; j++) {
    for(i = 1; i < imax+1; i++) {
      fprintf(fp, "%f\n", P[i][j] );
    }
  }

  if( fclose(fp) )
  {
    char szBuff[80];
    sprintf( szBuff, "Failed to close %s", szFileName );
    ERROR( szBuff );
  }
}


void write_vtkHeader( FILE *fp, int imax, int jmax, 
                      double dx, double dy) {
  if( fp == NULL )		       
  {
    char szBuff[80];
    sprintf( szBuff, "Null pointer in write_vtkHeader" );
    ERROR( szBuff );
    return;
  }

  fprintf(fp,"# vtk DataFile Version 2.0\n");
  fprintf(fp,"generated by CFD-lab course output (written by Tobias Neckel) \n");
  fprintf(fp,"ASCII\n");
  fprintf(fp,"\n");	
  fprintf(fp,"DATASET STRUCTURED_GRID\n");
  fprintf(fp,"DIMENSIONS  %i %i 1 \n", imax+1, jmax+1);
  fprintf(fp,"POINTS %i float\n", (imax+1)*(jmax+1) );
  fprintf(fp,"\n");
}


void write_vtkPointCoordinates( FILE *fp, int imax, int jmax,
                      double dx, double dy) {
  double originX = 0.0;
  double originY = 0.0;

  int i = 0;
  int j = 0;

  for(j = 0; j < jmax+1; j++) {
    for(i = 0; i < imax+1; i++) {
      fprintf(fp, "%f %f 0\n", originX+(i*dx), originY+(j*dy) );
    }
  }
}


