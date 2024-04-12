import os
import tempfile

def read_fasta(file_path):
    sequences = {}        #create a empty dictionary to store sequence
    with open(file_path, 'r') as file:   #open a file by 'read' method
        current_id = None
        current_sequence = ''
        for line in file:
            line = line.strip()          #make sure there is no any space 
            if line.startswith('>'):     #the fasta format start with >
                if current_id:           #make sure the sequence have id
                    sequences[current_id] = current_sequence
                current_id = line[1:]    #make sequence id  
                current_sequence = ''
            else:
                current_sequence += line #add all character into current_sequence
        if current_id:
            sequences[current_id] = current_sequence #make the key of sequence
    return sequences

def reverse_complement(sequence):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    reverse_sequence = sequence[::-1] #Use slicing to reverse order of the sequence from start to end
    reverse_complement_sequence = ''.join([complement[base] for base in reverse_sequence]) #make complementary sequence
    return reverse_complement_sequence

def main():
    file_path = input("Please input the path of the fasta file") #Ask the path of the fasta file
    sequences = read_fasta(file_path)   #use read_fasta function to read the file

    # Display sequence lengths
    print("Sequence lengths:")
    for identifier, sequence in sequences.items(): #identifier to store keys, sequence to store value
        print(f"{identifier}: {len(sequence)}") #print sequence ID and the length
    # Ask the user if they want to convert the file to reverse sequence
    convert_to_reverse = input("Do you want to convert the file to reverse sequence? (yes/no): ")
    if convert_to_reverse.lower() in ['yes', 'y']:
        # Create a temporary file and write the reverse complement sequences to it
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as temp_file:#use write mode, and don't delete file after close
            for identifier, sequence in sequences.items():
                reverse_complement_sequence = reverse_complement(sequence) #Use reverse function to make reverse sequence
                temp_file.write(f">{identifier}_reverse_complement\n") #write the ID of reverse sequence
                temp_file.write(f"{reverse_complement_sequence}\n")    #write the sequence of reverse sequence
                print(f"original sequence  ({identifier}): {sequence}")
                print(f"reverse sequence  ({identifier}): {reverse_complement_sequence}")

        # Request the user to input the location and the name for the new file.
        output_file_location = input("Please enter the Path for the new file: ")
        output_file_name = input("Please enter the name for the new file: ")
        output_file_path = os.path.join(output_file_location, output_file_name)

    # Move the temporary file to the location specified by the user.
        os.rename(temp_file.name, output_file_path)
        print("Reverse complement sequences have been successfully written to the new fasta file.")
    elif convert_to_reverse.lower() == 'no':
        print("No conversion performed. Exiting program.")
    else:
        print("Invalid input. Exiting program.")

if __name__ == "__main__":
    main()

