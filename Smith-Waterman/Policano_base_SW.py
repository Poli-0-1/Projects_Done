##############################################################

## Library needed to run the code:
import numpy as np
import pandas as pd             
from enum import IntEnum


## Cycle to run the code multiple times
while True:
        
    print("\nPlease insert the sequences:")         # Input seq 1 and seq 2

    Sequence_1 = input("\nSeq 1 : \n").strip()
    Sequence_2 = input("\nSeq 2 : \n").strip()


    sequen_1 = " " + Sequence_1.upper()       # Add space at the beginning of the sequence to correctly build the dataframe and uppercase for the correct format of the sequences
    sequen_2 = " " + Sequence_2.upper()
    
    def classic_style(seq_1, seq_2):            # Main function to align sequences 
        
        seq_A = sequen_1
        seq_B = sequen_2
        
        seq_A = [char for char in seq_A]
        seq_B = [char for char in seq_B]

        df_seq = pd.DataFrame(data = np.zeros((len(seq_B), len(seq_A))),
                            columns = seq_A,
                            index = seq_B)

        df_traceback = pd.DataFrame(data = np.zeros((len(seq_B), len(seq_A))),
                                    columns = seq_A,
                                    index = seq_B)


        Answer_scores = input("Would you like to change the scoring system? [Y/N] \n").upper()  # Ask the user if they want to change the scoring system
        if Answer_scores in ["Y", "y", "Yes"]:                # If the user wants to change the scoring system
            
             class SCORE(IntEnum):       # Scoring system for the alignment
                match = input("Please insert the score for a match: \n")
                mismatch = input("Please insert the score for a mismatch: \n")
                gap = input("Please insert the score for a gap: \n")
        else:                                   # If the user does not want to change the scoring system, use the 

            class SCORE(IntEnum):       # Scoring system for the alignment
                match = +1
                mismatch = -1
                gap = -1


        class Movement():                                   # Movement symbols for the traceback matrix
            mov_left = "\u2190"     # ← = Left-ward arrow
            mov_up = "\u2191"       # ↑ = Up-ward arrow
            mov_diag = "\\"      # ↖ = Diagonal arrow
            mov_stop = "\u00d7"     # ꘎ = Stop symbol

        maximum_score = 0                   # Keep track of the maximum score found during the alignment
        max_score_position = None           # Keep track of the position of the maximum score

        n_row, n_cols = df_seq.shape

        for i in range(1, n_row):                   # Cycle through the dataframe composed of the two sequence s
            for j in range(1, n_cols):
                if df_seq.columns[j] == df_seq.index[i]:
                    match_score = SCORE.match
                else:
                    match_score = SCORE.mismatch

                diag_score = df_seq.iloc[i - 1, j - 1] + match_score        # Assign score to the matrix 
                vert_score = df_seq.iloc[i - 1, j] + SCORE.gap
                later_score = df_seq.iloc[i, j - 1] + SCORE.gap

                df_seq.iloc[i, j] = max(0, diag_score, vert_score, later_score)

                if df_seq.iloc[i, j] >= maximum_score:                      # Update pf the maximum score and its position
                    maximum_score = df_seq.iloc[i, j]
                    max_score_position = (i, j)

                if df_seq.iloc[i, j] == 0:                                  # Assign the symbols to the traceback matrix      
                    df_traceback.iloc[i, j] = Movement.mov_stop
                elif df_seq.iloc[i, j] == diag_score:
                    df_traceback.iloc[i, j] = Movement.mov_diag
                elif df_seq.iloc[i, j] == later_score:
                    df_traceback.iloc[i, j] = Movement.mov_left
                elif df_seq.iloc[i, j] == vert_score:
                    df_traceback.iloc[i, j] = Movement.mov_up


        
        
        print("\nScore Matrix:\n")                  # Print the relevant information such as : "Score Matrix", "Traceback Matrix" 
        print(df_seq)
        
        
        print("\n---------------------------------------------------------------------------------------------------")    # Separator for better readability
        print("---------------------------------------------------------------------------------------------------\n")    # Separator for better readability
        
        
        print("\nTraceback Matrix:\n")
        print(df_traceback)

        print("\n---------------------------------------------------------------------------------------------------")    # Separator for better readability
        print("---------------------------------------------------------------------------------------------------\n")    # Separator for better readability
        
        
        seq_A_align = []                        # Section dedicated t the corect building of the alignement, had to figure out how to reverse the sequences and correctly adjust the construction logic
        seq_B_align = []

        i_max, j_max = max_score_position

        while df_traceback.iloc[i_max, j_max] != 0 and df_traceback.iloc[i_max, j_max] != Movement.mov_stop:        # Construction Logic of the alignment, using the traceback matrix to find the path taken during the alignment process
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

        seq_A_final = "".join(reversed(seq_A_align))        # Join the elements of the alignment lists and reverse them to get the correct order
        seq_B_final = "".join(reversed(seq_B_align))

        print("\nFinal Alignment:\n")         # Print the final alignment of the sequences
            
        print(f"\t{seq_A_final}")
        print(f"\t{seq_B_final}\n")
        print(f"\nHighest Score: \t {maximum_score}")
        
        return seq_A_final, seq_B_final

    classic_style(sequen_1, sequen_2)
                                                                # Ask the user if they want to try again        
    print("\nWould you like to try again?" + "\t     [Y/N]")
    
    answer = input().upper()
    if answer != "Y":
        break