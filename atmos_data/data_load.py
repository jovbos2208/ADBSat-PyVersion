import numpy as np
import pymsis as ps
import pandas as pd

column_names = ["MASS_DENSITY", "N2", "O2", "O", "HE", "H", "AR", "N", "AO", "NO", "TEMP"]

dates = np.arange(np.datetime64("2006-01-01T00:00"), np.datetime64("2007-01-01T00:00"), np.timedelta64(480, "m"))

lats = np.linspace(-90,90,13)
longs = np.linspace(-180,180,13)
heights = np.linspace(100,400,31)


data = ps.calculate(dates,longs,lats,heights)
for k,h in enumerate(heights):
    df = pd.DataFrame(columns=column_names)
    for i in range(13):
        for j in range(13):
            new_df = pd.DataFrame(data[:,i,j,k,:],columns=column_names)
            df = pd.concat([df,new_df], ignore_index=True)
    
    df.to_csv(f"atmos_data/database_{int(h):03d}km.csv",index=False)

    print(f'Dowloading data for alt: {int(h):03d}km done!')