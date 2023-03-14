import pandas as pd

#Export data 
def exportDataframe(downhole_data, wpermits, procdir):
    todayDate = datetime.date.today()
    todayDateStr = str(todayDate)
    
    downholeOutfile = procdir+r'\Downhole_BedrockPicks_'+todayDateStr+'.csv'
    wPermitOutFile = procdir+r'\wPermits_BedrockPicks_'+todayDateStr+'.csv'
    
    print('Exported '+downholeOutfile)
    print('Exported '+wPermitOutFile)
    
    downhole_data.to_csv(downholeOutfile, index_label="ID")
    wpermits.to_csv(wPermitOutFile, index_label="ID")