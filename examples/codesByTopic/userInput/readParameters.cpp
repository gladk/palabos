/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2017 FlowKit Sarl
 * Route d'Oron 2
 * 1010 Lausanne, Switzerland
 * E-mail contact: contact@flowkit.com
 *
 * The most recent release of Palabos can be downloaded at 
 * <http://www.palabos.org/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "palabos2D.h"
#include "palabos2D.hh"

using namespace plb;
using namespace std;
typedef double T;

int main(int argc, char* argv[]) {
    plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmp/");

    try {

    // 1. Read a parameter from the command line.
    int value = 27;
    global::argv(1).readNoThrow(value);
    pcout << value << endl;

    // 2. Read parameters from an XML file.
    string fName("demo.xml");

    XMLreader document(fName);

    string text;
    document["MyApp"]["Messages"]["Welcome"].read(text);
    pcout << "TEXT:" << text << endl;

    int number;
    bool boolVal;
    T someFloat;
    std::vector<int> numbers;

    document["MyApp"]["Value"]["OneNumber"].read(number);
    pcout << "ONE NUMBER:" << number << endl;

    document["MyApp"]["OneBool"].read(boolVal);
    pcout << "ONE BOOL:" << boolVal << endl;

    document["MyApp"]["OneFloat"].read(someFloat);
    pcout << "ONE FLOAT:" << someFloat << endl;

    document["MyApp"]["TwoNumbers"].read(numbers);
    pcout << "TWO NUMBERS:" << numbers[0] << " " << numbers[1] << endl;

    double nestedVal;
    XMLreaderProxy nested = document["MyApp"]["Nested"];
    for (; nested.isValid(); nested = nested.iterId()) {
        nested["NestedData"].read(nestedVal);
        pcout << "Nested value at id " << nested.getId() << ": " << nestedVal << endl;
    }


    // 3. Write data into an XML file.
    //
    XMLwriter demo2;
    XMLwriter& myApp = demo2["MyApp"];

    XMLwriter& messages = myApp["Messages"];
    messages.setString("Text of message");
    messages["Welcome"].setString("Welcome to MyApp");
    messages["FareWell"].setString("Good bye");
    XMLwriter& messages2 = myApp["Messages"][2];
    messages2["Welcome"].setString("Welcome again");

    vector<int> twoNumbers;
    twoNumbers.push_back(12);
    twoNumbers.push_back(24);
    myApp["Windows"]["Window"];
    myApp["OneNumber"][3].set(12);
    myApp["OneNumber"][2].set(2);
    myApp["OneNumber"][20].set(20);
    myApp["TwoNumbers"].set(twoNumbers);
    myApp["OneFloat"].set(12.34);

    demo2.print("demo2.xml");
    

    }

    catch (PlbIOException& exception) {
        pcout << exception.what() << endl;
        return -1;
    }
}
