from __future__ import division
from psychopy import visual,data,event,core,gui
import time
import numpy as np
from numpy import random
import os
import math
import pandas

#define amounts of everything
Nstim=2
Nresp=2
Ntrials=50

#Initialize data
response=0

print('initialization_ok')

#make data file
info= {"ppnr": 0}

# Data file
already_exists = True
while already_exists:
    myDlg = gui.DlgFromDict(dictionary = info, title = "test")
    directory_to_write_to = "/Users/Pieter/Documents/psychopy_exercices" + "/" + str(info["ppnr"])
    if not os.path.isdir(directory_to_write_to):
        os.mkdir(directory_to_write_to)
    file_name = directory_to_write_to + "/test" + str(info["ppnr"]) 
    if not os.path.isfile(file_name+".tsv"):
        already_exists = False
    else:
        myDlg2 = gui.Dlg(title = "Error")
        myDlg2.addText("Try another participant number")
        myDlg2.show()
print("OK, let's get started!")

thisExp = data.ExperimentHandler(dataFileName = file_name, extraInfo = info)

#clocks
my_clock= core.Clock()
sleepclock=core.Clock()

#make window
window=visual.Window(fullscr=True, monitor='testMonitor', color=[0,0,0],colorSpace='rgb') #size=(800,600),

print('window_ok')



#elements of trials
t1 = [0.1]*Ntrials
t2= [0.5]*Ntrials
t3=[1]*Ntrials
sleeptime=np.column_stack((t1, t2, t3))
print(sleeptime)

#stimuli for later
fixation = visual.TextStim(window,text=("+"))

vert_grating=visual.GratingStim(win=window,tex='sin', mask='raisedCos',ori=45, pos=(0,0), size= (7,7), units='deg', interpolate=False)
gratings=[vert_grating]

print('graphical_elements_ok')

#randomize stimulus presentations (equal amounts per part)
stim=np.zeros(Ntrials)

s=np.repeat(np.arange(Ngratings),(Ntrials/Ngratings))
s=np.random.permutation(s)
stim=s

print('stim_ok')

# Within-subjects design
TrialList=[]
TrialList.append( {"sleeptimes": sleeptime})
print(TrialList)
trials = data.TrialHandler(trialList = TrialList, nReps=1, method = "sequential")
thisExp.addLoop(trials)

# welcome screen
welcome=visual.TextStim(window,text=("Welkom"))
welcome.draw()
window.flip()
event.waitKeys(keyList = "space")

sleepclock.reset()
#trial loop   
for trial in range(Ntrials):
    
    trials.addData("Tr", trial)
    
    t1_on=sleepclock.getTime()
    time.sleep(sleeptime[0,trial])
    t1_off=sleepclock.getTime()
    t1_diff=t1_off-t1_on
    
    t2_on=sleepclock.getTime()
    time.sleep(sleeptime[1,trial])
    t2_off=sleepclock.getTime()
    t2_diff=t2_off-t2_on
    
    t3_on=sleepclock.getTime()
    time.sleep(sleeptime[0,trial])
    t3_off=sleepclock.getTime()
    t3_diff=t3_off-t3_on
    
    trials.addData("t1_on", t1_on)
    trials.addData("t1_off", t1_off)
    trials.addData("t1_diff", t1_diff)
    
    trials.addData("t2_on", t2_on)
    trials.addData("t2_off", t2_off)
    trials.addData("t2_diff", t2_diff)
    
    trials.addData("t3_on", t3_on)
    trials.addData("t3_off", t3_off)
    trials.addData("t3_diff", t3_diff)

#say goodbye to the participant
goodbye         = visual.TextStim(window,text=( "Dada")),wrapWidth=50, units="deg",alignHoriz='center', alignVert='bottom')
goodbye.draw()
window.flip()
event.waitKeys(keyList = "space")

thisExp.saveAsWideText(file_name, appendFile=False)
thisExp.abort()
window.close()
core.quit()    