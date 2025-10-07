####################################################################################

### For semplicity reasons the comments of the modified version of the code will only be present in the effectively modified section
## If you want to examine the code with the comments, please refer to the original file: final_base_SW.py

import numpy as np
import pandas as pd
from enum import IntEnum


while True:
    
    print("\nPlease insert the sequences:")

    Sequence_1 = input("\nSeq 1 : \n").strip()
    Sequence_2 = input("\nSeq 2 : \n").strip()


    sequen_1 = " " + Sequence_1.upper()
    sequen_2 = " " + Sequence_2.upper()

    def modified_style(seq_1, seq_2):
        seq_A = [char for char in seq_1]
        seq_B = [char for char in seq_2]

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
                

        class Movement():
            mov_left = "\u2190"     # ← = Left-ward arrow
            mov_up = "\u2191"       # ↑ = Up-ward arrow
            mov_diag = "\\"         # ↖ = Diagonal arrow
            mov_stop = "\u00d7"     # ꘎ = Stop symbol

        maximum_score = 0
        max_score_position = None

        n_row, n_cols = df_seq.shape

        for i in range(1, n_row):
            for j in range(1, n_cols):
                if df_seq.columns[j] == df_seq.index[i]:
                    match_score = SCORE.match
                else:
                    match_score = SCORE.mismatch

                diag_score = df_seq.iloc[i - 1, j - 1] + match_score
                vert_score = df_seq.iloc[i - 1, j] + SCORE.gap
                later_score = df_seq.iloc[i, j - 1] + SCORE.gap

                df_seq.iloc[i, j] = max(0, diag_score, vert_score, later_score)

                if df_seq.iloc[i, j] >= maximum_score:
                    maximum_score = df_seq.iloc[i, j]
                    max_score_position = (i, j)

                if df_seq.iloc[i, j] == 0:
                    df_traceback.iloc[i, j] = Movement.mov_stop
                elif df_seq.iloc[i, j] == diag_score:
                    df_traceback.iloc[i, j] = Movement.mov_diag
                elif df_seq.iloc[i, j] == later_score:
                    df_traceback.iloc[i, j] = Movement.mov_left
                elif df_seq.iloc[i, j] == vert_score:
                    df_traceback.iloc[i, j] = Movement.mov_up


        
        
        print("\nScore Matrix:\n")
        print(df_seq)
        
        print("\n---------------------------------------------------------------------------------------------------")    # Separator for better readability
        print("---------------------------------------------------------------------------------------------------\n")    # Separator for better readability
        
        
        print("\nTraceback Matrix:\n")
        print(df_traceback)

        print("\n---------------------------------------------------------------------------------------------------")    # Separator for better readability
        print("---------------------------------------------------------------------------------------------------\n")    # Separator for better readability
        

        
        
        seq_A_align = []
        seq_B_align = []
        
        i_max, j_max = max_score_position
        
        
        
        one_gap_inside = []  # List to store the sequences with at least one three consecutive matches
        
        
        threshold = (80 * maximum_score)/100       # Threshold for the alignement score, those elements that do not reach the minimum value will not be printed in the final output
        
        
        
        while df_traceback.iloc[i_max, j_max] != 0 and df_traceback.iloc[i_max, j_max] != Movement.mov_stop:
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

        seq_A_final = "".join(reversed(seq_A_align))
        seq_B_final = "".join(reversed(seq_B_align))
            


        def One_Only_Gap(seq_A_final, seq_B_final):        # Function dedicated to find the three consecutives mathces in the sequences 
            
            Inside_Gap = False
            
            One_Gap = 0                                          #Implement a cycle in order to count the times we have a match between the two sequences and the 3 consecutive matches                    
            Longest_Gap = 0
            Gap_position = 0
            
            for a, b in zip(seq_A_final, seq_B_final):                  # Zip the two sequences together to compare them
                if a == "-" or b == "-":
                    
                    if not Inside_Gap:                          # If we find a gap, we set the Is_a_Gap variable to True
                        Inside_Gap = True
                        One_Gap += 1                             # add "points" to the variables or reset them 
                        Gap_position = 1
                    
                    else:
                        Gap_position += 1


                        
                        
                else :
                    if Inside_Gap:  
                        Longest_Gap = max(Longest_Gap, Gap_position)  # If we find a match, we set the Is_a_Gap variable to False and reset the Gap_position
                        
                        Inside_Gap = False      # Reset the Inside_Gap variable to False if we find a match
                        Gap_position = 0

                    
                                
            if Inside_Gap:
                Longest_Gap = max(Longest_Gap, Gap_position)    # check the last gap if the loop ends with a gap 
                    
                
            return Longest_Gap, One_Gap
        
        
        
        Longest_Gap, One_Gap = One_Only_Gap(seq_A_final, seq_B_final)  # Apply the function to the sequences
        
        max_score_alignment = df_seq.iloc[max_score_position]
                    

        
        if One_Gap == 1 and max_score_alignment >= threshold:  # If the sequences have at least one gap and the score is above the threshold, we append them to the list
            one_gap_inside.append((seq_A_final, seq_B_final, max_score_alignment, Longest_Gap))
                        
        
        unique_align = list(set(one_gap_inside))                                          # Find the unique alignement of the sequences that have at least one three consecutive matches
        sorted_final_align = sorted(unique_align, key = lambda x : x[3], reverse = True)    # Sort them based on the score, in descending order
     
     
        if sorted_final_align:
            highest_gap_sequences = sorted_final_align[0]  # Get the first element of the sorted list, which is the one with the highest score

            print(f"\nHighest Gap Scoring Alignment Sequences:\n{highest_gap_sequences[0]}\n{highest_gap_sequences[1]}\n")
            
            print(f"\nHighest Scoring: {highest_gap_sequences[2]}\n")  # Print the score of the alignment
            
            print(f"\nLongest Gap in the sequences:\n {highest_gap_sequences[3]}\n")  # Print the longest gap in the sequences
        
            print("\n--------------------------------------------\n")
        
            if len(sorted_final_align) > 1:
            
                sorted_final_align.sort(key=lambda x: x[2], reverse = False)  # Sort the alignments by score in ascending order
                for h, (s1, s2, score, longest_gap) in enumerate(sorted_final_align[1:],1):
                    
                    print(f"{h}. Alignment with Score: {score}, Longest Gap: {longest_gap}")
                    print(s1)
                    print(s2)
                    
                    print("\n--------------------------------------------\n")

            
    modified_style(sequen_1, sequen_2)
    
    print("\nWould you like to try again? " + "\t   \U0001f408  [Y/N]")
    answer = input().upper()
    if answer != "Y":
        break