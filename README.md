
################################################################################
# JointAngles - README
################################################################################
# 1. info
################################################################################
#
# If you use this software, please cite it as below:
#
# authors:
# - family-names: "Campo"
#  given-names: "Adriaan"
#  orcid: "https://orcid.org/0000-0003-0219-6054"
# title: "Extracting joint angles from MoCap Data: a Matlab package."
# version: 1.0.0-alpha
#
# date-released: 2023-04-04
#
################################################################################
# 2. readme
################################################################################
#
# This code will calculate joint angles from the Animation Marker Set, as used
# by Qualisys. The code is a MATLAB class, and takes solely the path to a MoCap
# data file as an input, in *.mat format.
#
# To calculate joint angles run the code as follows in the MATLAB command
# window:
#
# data = ReadMocapData('full file path')
# data = data.PreProcessData;
# data = data.ProcessData;
#
# data now contains time series of all the joint angles, for the left and right
# side of the body, organised per joint, and per measured angle:
# data.Participants.ProcessedData.L.JointAngleData.[Joint].[Angle], for the left
# side and
# data.Participants.ProcessedData.R.JointAngleData.[Joint].[Angle], for the
# right side.
#
# Evaluated joints/angles are:
# Elbow: flexion (FE), pronation (PS)
# Wrist: adduction (AA), flexion (FE)
# Shoulder: adduction (AA), elevation (E), pronation (PS)
# Ankle: dorsiflexion (DP), elevation (E), rotation (R)
# Knee: flexion (FE)
# Hip: adduction (AA), flextion (FE)
# Neck: alpha (AA), beta (BB), gamma (GG)
# Spine: alpha (AA), beta (BB), gamma (GG)
#
# for example, a column with all the angles of right elbow flexion would be
# generated as:
#
# FID = data.Participants.ProcessedData.R.JointAngleData.Elbow.FE.Angle
#
################################################################################
# 3. contact
################################################################################
#
# This code is under constant development. For questions, comments and concerns,
# you can contact me at: adriaan.campo@ugent.be
#
################################################################################
