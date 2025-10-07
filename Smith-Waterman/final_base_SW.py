import numpy as np
import pandas as pd
from enum import IntEnum


while True:                                     # Use of a cycle to repeat the process in case of changes from the use in the input sequence
        
    print("\nPlease insert the sequences:")

    Sequence_1 = input("\nSeq 1 : \n").strip()          # Reading and correct formatting of the sequences
    Sequence_2 = input("\nSeq 2 : \n").strip()


    sequen_1 = " " + Sequence_1.upper()                         
    sequen_2 = " " + Sequence_2.upper()
    
    def classic_style(seq_1, seq_2):                    # Function for the SW algorithm
        
        seq_A = sequen_1
        seq_B = sequen_2
        
        seq_A = [char for char in seq_A]            # Reading of the sequenc's elements one by one 
        seq_B = [char for char in seq_B]

        df_seq = pd.DataFrame(data = np.zeros((len(seq_B), len(seq_A))),        # Construction of dataframe
                            columns = seq_A,
                            index = seq_B)

        df_traceback = pd.DataFrame(data = np.zeros((len(seq_B), len(seq_A))),  # Construction of traceback dataframe
                                    columns = seq_A,
                                    index = seq_B)

        class SCORE(IntEnum):        # Scores settings
            match = +1
            mismatch = -1
            gap = -1

        class Movement():
            mov_left = "\u2190"     # ← = Left-ward arrow
            mov_up = "\u2191"       # ↑ = Up-ward arrow
            mov_diag = "\\"         
            mov_stop = "\u00d7"     # ꘎ = Stop symbol

        maximum_score = 0               # Max score and its position to later use in the traceback
        max_score_position = None

        n_row, n_cols = df_seq.shape

        for i in range(1, n_row):                  # Assign values depending on the elements considered between the two sequences
            for j in range(1, n_cols):
                if df_seq.columns[j] == df_seq.index[i]:
                    match_score = SCORE.match
                else:
                    match_score = SCORE.mismatch

                diag_score = df_seq.iloc[i - 1, j - 1] + match_score        # Assign the scores to later use them in the function
                vert_score = df_seq.iloc[i - 1, j] + SCORE.gap
                later_score = df_seq.iloc[i, j - 1] + SCORE.gap

                df_seq.iloc[i, j] = max(0, diag_score, vert_score, later_score)     # Find the maximum value

                if df_seq.iloc[i, j] >= maximum_score:          # If found, substitute the max score
                    maximum_score = df_seq.iloc[i, j]
                    max_score_position = (i, j)

                if df_seq.iloc[i, j] == 0:                          # Apply the obtained value in the traceback dataset for easier backward movement
                    df_traceback.iloc[i, j] = Movement.mov_stop
                elif df_seq.iloc[i, j] == diag_score:
                    df_traceback.iloc[i, j] = Movement.mov_diag
                elif df_seq.iloc[i, j] == later_score:
                    df_traceback.iloc[i, j] = Movement.mov_left
                elif df_seq.iloc[i, j] == vert_score:
                    df_traceback.iloc[i, j] = Movement.mov_up


        
        
        print("\nScore Matrix:")            # Print to the user both Matrices (Normal and Traceback)
        print(df_seq)
        print("\nTraceback Matrix:")
        print(df_traceback)

        
        seq_A_align = []
        seq_B_align = []

        i_max, j_max = max_score_position       # Identify the i, j (position) for the maximum score

        while df_traceback.iloc[i_max, j_max] != 0 and df_traceback.iloc[i_max, j_max] != Movement.mov_stop:        # Compute the traceback while assigning and using points
            if df_traceback.iloc[i_max, j_max] == Movement.mov_diag:
                seq_A_align.append(seq_A[j_max])
                seq_B_align.append(seq_B[i_max])
                i_max -= 1
                j_max -= 1
            elif df_traceback.iloc[i_max, j_max] == Movement.mov_left:
                seq_A_align.append(seq_A[j_max])
                seq_B_align.append("-")
                j_max -= 1
            elif df_traceback.iloc[i_max, j_max] == Movement.mov_up:
                seq_A_align.append("-")
                seq_B_align.append(seq_B[i_max])
                i_max -= 1

        seq_A_final = "".join(reversed(seq_A_align))        # Compute the final sequences A and B (which we have to reverse as they are backward in the traceback matrix)
        seq_B_final = "".join(reversed(seq_B_align))

        print("\nFinal Alignment:")         # Print the final alignment
            
        print(seq_A_final)                  # Print sequence A and B
        print(seq_B_final)
        return seq_A_final, seq_B_final

    classic_style(sequen_1, sequen_2)

    print("Would you like to try again?  [Y/N]")        # Cycle to try again
    answer = input().upper()
    if answer != "Y":
        break