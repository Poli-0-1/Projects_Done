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

        class SCORE(IntEnum):
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


        
        
        print("Score Matrix:")
        print(df_seq)
        print("\nTraceback Matrix:")
        print(df_traceback)

        
        seq_A_align = []
        seq_B_align = []
        
        i_max, j_max = max_score_position
        
        
        
        three_of_a_match = []

        threshold = (60 * maximum_score)/100
        
        
        
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

        print("\nFinal Alignment:")
            
        print(seq_A_final)
        print(seq_B_final)

        def three_consecutive_matches(seq_A_final, seq_B_final):
            
            match_count = 0                                          #Use an iterative cycle to find the times we have a match
            occurences = 0                                           # and at how many match we count it as occurance for the filtering
            
            for a,b in zip(seq_A_final, seq_B_final):
                if a == b and a != "-" and b != "-":
                    match_count += 1   
                    if match_count == 3:                               #if the matches chain is disrupted return the count to 0 and start again
                        occurences += 1
                        match_count = 0
                else :  
                    match_count = 0
            return occurences
        
        
        def longest_streak(seq_A_final, seq_B_final):
            streak = 0
            longest_streak = 0

            for a,b in zip(seq_A_final, seq_B_final):                   #same logic for the Three_of_a_batch function but with slight changes
                if a == b and a != "-" and b != "-":
                    streak += 1
                    if streak > longest_streak:
                        longest_streak = streak        
                else :
                    streak = 0
                    
            return longest_streak
        
        
        for i in range(1, n_row):
            for j in range(1, n_cols):
                if df_seq.iloc[i, j] >= threshold:
                    
                    if three_consecutive_matches(seq_A_final, seq_B_final) >= 1:
                        three_of_a_match.append((seq_A_final, seq_B_final, df_seq.iloc[i, j]))
                
        
        unique_align = list(set(three_of_a_match))
        sorted_final_align = sorted(unique_align, key = lambda x : longest_streak(x[0], x[1]), reverse = True)
        
        three_of_a_match_value = three_consecutive_matches(seq_A_final, seq_B_final)
        longest_streak_value = longest_streak(seq_A_final, seq_B_final)
        
        
        print(sorted_final_align)
        print(f"\n 3 consecutive matches count:\n {three_of_a_match_value}")
        print(f"\n Longest streak:\n {longest_streak_value}")
        
    modified_style(sequen_1, sequen_2)
    
    print("Would you like to try again?  [Y/N]")
    answer = input().upper()
    if answer != "Y":
        break