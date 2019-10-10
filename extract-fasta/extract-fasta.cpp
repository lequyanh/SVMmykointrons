#include <iostream>
#include <unordered_map>
#include <sstream>
#include <assert.h>
#include <fstream>
#include <algorithm>
#include <string.h>

using namespace std;

const string NEGATIVE_STRAND = "-";
const string POSITIVE_STRAND = "+";

void print_usage() {
    cout << "Usage: extract-fasta [-i] FASTA" << endl;
    cout << "Extract sequences from a FASTA file." << endl;
    cout << endl;
    cout << "Use option -i to make the program ignore incomplete subsequences (they will not be outputed)." << endl;
    cout << endl;
    cout << "The input lines are expected to be formatted as follows:" << endl;
    cout << "SCAFFOLD STRAND START1 END1 START2 END2 ... STARTn ENDn" << endl;
    cout << endl;
    cout << "SCAFFOLD is the name of scaffold in the FASTA." << endl;
    cout << "STRAND   is the strand (" << POSITIVE_STRAND << " or " << NEGATIVE_STRAND << ")" << endl;
    cout << "STARTi   is the start position (inclusively) of i-th subsequence within the scaffold." << endl;
    cout << "ENDi     is the end position (inclusively) of i-th subsequence within the scaffold." << endl;
    cout << endl;
    cout << "At least one pair of start and end positions must be specified." << endl;
    cout << endl;
    cout << "Example input:" << endl;
    cout << "chromosome_5 + 1 30 43 71 100 109" << endl;
    cout << "chromosome_3 - 25 34 89 99" << endl;
    cout << "chromosome_1 + 35 45" << endl;
}

unordered_map<string, string> load_data(istream &input_stream) {
    unordered_map<string, string> data;
    string current_scaffold;
    stringstream seq_buffer;

    while (true) {
        input_stream >> ws; //eat all white spaces

        int current_char = input_stream.peek();

        if (current_char == EOF) {
            break;
        } else if (current_char == '>') {
            //check if there are some data for previous scaffold before moving to next one
            string sequence = seq_buffer.str();
            if (!sequence.empty()) {
                //current scaffold should not be empty
                assert(!current_scaffold.empty());
                //there must not be the same scaffold in the map already
                //that would mean that the same scaffold is defined twice
                assert(data[current_scaffold].empty());

                //store the scaffold
                data[current_scaffold] = sequence;

                //clear the buffer
                seq_buffer.str("");
            }

            //ignore '>' since it is not part of the scaffold name
            input_stream.ignore(1);

            //read the scaffold name
            input_stream >> current_scaffold;
        } else {
            //read sequence into the buffer
            string sequence;
            getline(input_stream, sequence);
            seq_buffer << sequence;
        }
    }

    //store the last remaining scaffold sequence
    string sequence = seq_buffer.str();
    if (!sequence.empty()) {
        assert(!current_scaffold.empty());
        assert(data[current_scaffold].empty());
        data[current_scaffold] = sequence;
    }

    return data;
}

void make_complementary(string &sequence) {
    for (int i = 0; i < sequence.size(); ++i) {
        switch (sequence[i]) {
            case 'A':
            case 'a':
                sequence[i] = 'T';
                break;
            case 'T':
            case 't':
                sequence[i] = 'A';
                break;
            case 'C':
            case 'c':
                sequence[i] = 'G';
                break;
            case 'G':
            case 'g':
                sequence[i] = 'C';
                break;
        }
    }
    reverse(begin(sequence), end(sequence));
}

int main(int argc, const char **argv) {
    bool ignore_incomplete = false;
    if (argc < 2 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
        print_usage();
        return 0;
    } else if (argc > 2 && strcmp(argv[1], "-i") == 0) {
        ignore_incomplete = true;
    }

    //turn off synchronization with C streams (speed up)
    ios_base::sync_with_stdio(false);

    ifstream input_file(argv[argc - 1]);

    //maps name of a scaffold to its sequence
    unordered_map<string, string> data = load_data(input_file);


    string line;

    //read whole line
    next_line: while (getline(cin, line)) {
        istringstream input_line(line);
        ostringstream output_sequence_stream;

        string scaffold, strand;
        int start, end;

        //read scaffold and strand
        input_line >> scaffold >> strand;

        //read individual start-end pairs
        while ((input_line >> start >> end)) {
            //start must be indexed from 1
            assert(start > 0);
            assert(strand == NEGATIVE_STRAND || strand == POSITIVE_STRAND);
            //retrieve a sequence for the scaffold
            string &scaffold_sequence = data[scaffold];

            //the sequence must be sufficiently long
            if (scaffold_sequence.size() < end) {
                if (ignore_incomplete) goto next_line;
                else {
                    cerr << "Incomplete sequence. Use -i to skip such sequences." << endl;
                    return 1;
                }
            }

            //extract the desired subsequence (string is, unlike start positions, indexed from 0 => start-1)
            string extracted = scaffold_sequence.substr(start - 1, end - start + 1);
            output_sequence_stream << extracted;
        }

        //print all the subsequences as a single sequence
        string output_sequence = output_sequence_stream.str();

        if (strand == NEGATIVE_STRAND) {
            make_complementary(output_sequence);
        }

        transform(output_sequence.begin(), output_sequence.end(), output_sequence.begin(), ::toupper);

        cout << '>' << line << endl;
        cout << output_sequence << endl;        
    }

    return 0;
}
