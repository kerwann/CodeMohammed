################################################################ 

## utilisation des variables temporelles
import time

## utilisation des regex
import re

## utilisation des fonctions numeriques
import numpy as np

## 
import imageio

## 
import os

##
from math import * 

## 
import glob

## utilisation des fonctions de lecture dicom
import pydicom 

## 
from natsort import natsorted

## 
from pydicom.uid import generate_uid

##
from scipy import interpolate

## Dataframe
import pandas as pd

## Graphique
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.patches as patches
from mpl_toolkits import mplot3d

## Gere les arguments
from argparse import ArgumentParser

##Utiliser pour lister les dossiers
import os

## Tiff
from PIL import Image

##
import cv2

##
import random as rng

#################################################################
## Function #####################################################
################################################################# 

def getExternalRoiNumber (RTStruct) : 
    """Get the ROI number for the External Structure

    Parameters:
    RT struct Dicom (pydicom): RT struct Dicom File
    
    Returns:
    int[]: Roi number 

    """
    roiList = len(RTStruct.StructureSetROISequence)
    ROInumber = -1
    for i in range ((roiList)) : 
        if RTStruct.StructureSetROISequence[i].ROIName=="External" :
            ROInumber = i
    return ROInumber

def getRoiPointCoordinates (RTstruct, Roinumber) : 
    """Get ROI point coordinates 

    Parameters:
    RT struct Dicom (pydicom): RT struct Dicom File
    Roinumber (int) : number of ROI
    
    Returns:
    List[]: List of ANp.Array containning points coordinates  

    """
    RTstructLength = len(RTstruct.ROIContourSequence[(Roinumber)].ContourSequence) 
 
    pointsCoordinatesXYZ = []
    z_coordinates = []
     
    for i in range(RTstructLength) :
        temp = [RTstruct.ROIContourSequence[Roinumber].ContourSequence[i].ContourData]
        temp = np.array(temp, dtype=float)
        temp = np.reshape(temp,(-1,3))
        
        # print(temp[0][2])
        if not temp[0][2] in z_coordinates : 
            z_coordinates.append(temp[0][2])
            pointsCoordinatesXYZ.append(temp)    
        else :
            index = z_coordinates.index(temp[0][2])
            pointsCoordinatesXYZ[index] = np.concatenate((pointsCoordinatesXYZ[index],temp),axis = 0)   
        
        if i==0 :
            minList = np.amin(temp, axis=0)
            maxList = np.amax(temp, axis=0)
            
        minListTemp = np.amin(temp, axis=0)
        maxListTemp = np.amax(temp, axis=0)
        for i in range (0,minListTemp.shape[0]) :
            if minListTemp[i] < minList[i] :
                minList[i] = minListTemp[i]
            if maxListTemp[i] > maxList[i] :
                maxList[i] = maxListTemp[i] 
                   
    return pointsCoordinatesXYZ,minList,maxList

def saveCoordinates (PointsCoordinates, directory) :
    """Save point coordinates 

    Parameters:
    PointCoordinates (List of np.array) : List of np.array by slice
    directory (str) : directory where to save file 
    
    """
    fileName = directory + "coordinates.txt"
    for l in range (0,len(PointsCoordinates)) :
        for row in PointsCoordinates[l] :
            fTxt = open(fileName, 'a')
            fTxt.write(str(row[0])+"\t"+str(row[1])+"\t"+str(row[2])+"\n")
            fTxt.close    

def getKey(item):
    return item[0]

def getDCMname(tupleList):
    """Return the filename."""
    orderNameList = []
    for i in tupleList:
        orderNameList.append(i[1])
    return orderNameList

def orderRTstructDir(path):
    """Order the dicom directory w.r.t the position

    Parameters:
    path (str): path to the dicom directory

    Returns:
    str[]: List containing the dicom filename ordered by the position 

   """
    position = [0]
    filelist = [file for file in os.listdir(os.path.join(path)) if file.endswith('.dcm') if not file.startswith(".")]
    return filelist

def orderDicomDir(path):
    """Order the dicom directory w.r.t the position

    Parameters:
    path (str): path to the dicom directory

    Returns:
    str[]: List containing the dicom filename ordered by the position 

   """
    position = []
    filelist = [file for file in os.listdir(os.path.join(path)) if file.endswith('.dcm') if not file.startswith(".")]
    for filename in filelist:
        ds = pydicom.dcmread(os.path.join(path,filename))
        if hasattr(ds, 'SliceLocation'):
            position.append((float(ds.SliceLocation),filename))
        else:
            position.append((float(ds.ImagePositionPatient[2]),filename))
    
    return getDCMname(sorted(position, key=getKey))

def getNumberOfSlice (DicomDir) : 
    
    dicom_list = orderDicomDir(DicomDir)
    
    return len(dicom_list)

def getCTinformation (DicomDir) :

    dicom_list = orderDicomDir(DicomDir)
    dsCT = pydicom.dcmread(os.path.join(DicomDir,dicom_list[0]))
    
    return dsCT
    
def getRTstructInformation (DicomDir) :

    dicom_list = orderRTstructDir(DicomDir)
    dsRTStruct = pydicom.dcmread(os.path.join(DicomDir,dicom_list[0]))
    
    return dsRTStruct 

def createExternalBasedMask(CTpath,RTstructPath) : 
    
    dsCT = getCTinformation(CTpath) 
    dsRTStruct = getRTstructInformation(RTstructPath)
    
    #### Recupere les coordonnees des points composants l'externe
    ROInumber = getExternalRoiNumber(dsRTStruct)
    coordinates, minList, maxList = getRoiPointCoordinates (dsRTStruct, ROInumber)
    
    #### Recupere les informations sur l'image
    
    ## Origine Dicom CT
    x_0 = dsCT.ImagePositionPatient[0]
    y_0 = dsCT.ImagePositionPatient[1]
    z_0 = dsCT.ImagePositionPatient[2]
    # print("Origine CT : (",x_0,";",y_0,";",z_0,")")
    
    ## Matrix Size
    nx = dsCT.Rows
    ny = dsCT.Columns
    nz = getNumberOfSlice (CTOriginial_path)
    # print("Taille de la matrice CT : (",nx,";",ny,";",nz,")")
    
    ## Resolution de la grille CT
    dx = dsCT.PixelSpacing[0]
    dy = dsCT.PixelSpacing[1]
    dz = dsCT.SliceThickness
    # print("Resolution CT : (",dx,";",dy,";",dz,")")
    
    ## Opposé de l'origine Dicom CT
    x_CTmax = dsCT.ImagePositionPatient[0] + nx * dx
    y_CTmax = dsCT.ImagePositionPatient[1] + ny * dy
    z_CTmax = dsCT.ImagePositionPatient[2] + nz * dz
    # print("Origine CT : (",x_0,";",y_0,";",z_0,")")
    
    #### Initialisation de la matrice pour la création du masque
    M = np.zeros((ny,nx,nz))
    # M = np.empty((nx,ny,nz), dtype=int)
    # M[:] = np.nan
    
    ## Bounding Box
    ## Axe X : Droite = négatif / Gauche = Positif
    ## Axe Y : Ant = négatif / Post = Positif
    ## Axe Z : Inf = négatif / Sup = Positif
  
    ## Limite dans le repere (x,y,z)    
    x_min = minList[0]
    y_min = minList[1]
    z_min = minList[2]
    # print(x_min,y_min,z_min)

    x_max = maxList[0]
    y_max = maxList[1]
    z_max = maxList[2]
    # print(x_max,y_max,z_max)

    PointList = []
    for y in np.arange(y_min,y_max,dy):
        for x in np.arange(x_min,x_max,dx):
            PointList.append([x,y])
    y_bins = int((y_max-y_min)/dy) + 1
    x_bins = int((x_max-x_min)/dx) + 1
    z_bins = len(coordinates)
    
    # print(y_bins, x_bins, z_bins)
    
    n_right = int(np.around((x_min-x_0)/dx))
    # n_left = int(np.around((x_CTmax-x_min)/dx))-x_bins
    n_left = nx - x_bins - n_right
    rightColumn_array = np.zeros((y_bins,n_right), dtype=bool)
    leftColumn_array = np.zeros((y_bins,n_left), dtype=bool)
 
    n_ant = int(np.around((y_min-y_0)/dy))
    # n_post = int(np.around((y_CTmax-y_min)/dy))-y_bins
    n_post = ny - y_bins - n_ant
    antRow_array = np.zeros((n_ant,nx), dtype=bool)
    postRow_array = np.zeros((n_post,nx), dtype=bool)

    n_inf = int(np.around((z_min-z_0)/dz))
    # n_sup = int(np.around((z_CTmax-z_min)/dz))-z_bins
    n_sup = nz - z_bins - n_inf
        
    supArray = np.zeros((ny,nx,n_sup), dtype=bool)
    infArray = np.zeros((ny,nx,n_inf), dtype=bool)

    #### Utilisation de Path de matplolib
    #### Attention, il faudra encore voir le cas ou il met 2 coupes par contours, voir au moment du remplissage de coordinates
        
    for r in range (0,len(coordinates)) :  
        points = coordinates[r]
        z = coordinates[r][0][2]
        points = np.delete(points,2,axis=1)    
        path_contour = mpath.Path(points)
        res = path_contour.contains_points(PointList)

        ## Mise en forme de la matrice 
        M = np.reshape(res,(y_bins,x_bins))
                
        ## Mise en forme en 512 x 512
        M = np.concatenate((rightColumn_array,M), axis = 1)
        M = np.concatenate((M,leftColumn_array), axis = 1)
        M = np.concatenate((antRow_array,M), axis = 0)
        M = np.concatenate((M,postRow_array), axis = 0)
        
        M = np.reshape(M,(M.shape[0],M.shape[1],1))
        
        if (r == 0) :
            Mask = M        
        else :
            Mask = np.concatenate((Mask,M), axis = 2)
    
    ## Concatener la matrice dans l'axe z-z_0
    Mask = np.concatenate((infArray,Mask), axis = 2)
    Mask = np.concatenate((Mask,supArray), axis = 2)
        
    return Mask
        
#################################################################
## Path #########################################################
################################################################# 

MotherPath = "..\Etude_mDD\CheeseModifInsert"
subdirList = ['Cheese_100kV','Cheese_120kV','Cheese_140kV','Cheese_80kV']
subdirList = ['Cheese_120kV']


MotherPath = "..\Etude_mDD\Patient"
subdirList = ['109127-Clinac']

for l in subdirList:

    Path = MotherPath+"\\"+l
    print(Path)  

    
    #################################################################
    ## DICOM PATH ###################################################
    #################################################################
    
    #### Image Original ###########
    CTOriginial_path = Path+'\Original'

    #### DICOM : RTstruct ###########
    RTstruct_path = Path+'\RTstruct'
    

    #################################################################
    ## Creation du masques ##########################################
    #################################################################
    
    print("Creation du masque : ")

    M = createExternalBasedMask(CTOriginial_path,RTstruct_path) 
    


    