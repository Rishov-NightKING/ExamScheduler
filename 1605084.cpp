#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>
#include <queue>
#include <set>
#include <cmath>
#include <tuple>
#include <unordered_set>
#include <algorithm>
#include <sstream>
#include <iterator>
#include <random>
#include <cstdlib>

using namespace std;

template <typename Iter, typename RandomGenerator>
Iter select_randomly(Iter start, Iter end, RandomGenerator &g)
{
    std::uniform_int_distribution<> dis(0, std::distance(start, end) - 1);
    std::advance(start, dis(g));
    return start;
}

template <typename Iter>
Iter select_randomly(Iter start, Iter end)
{
    static std::random_device rd;
    static std::mt19937 gen(rd());
    return select_randomly(start, end, gen);
}

int getRandomBetweenTwo(int a, int b)
{
    return (rand() > RAND_MAX / 2) ? a : b;
}

vector<int> kempeNodes;

class Graph
{

public:
    int vertices;
    unordered_set<int> *adjacencylist;
    vector<int> allVerticesVector;

    int *timeslots;
    bool *discovered;

    Graph(int V)
    {
        vertices = V;
        adjacencylist = new unordered_set<int>[V];
        timeslots = new int[V];
        discovered = new bool[V];
        
        for (int i = 0; i < V; i++)
        {
            allVerticesVector.push_back(i);
        }
    }

    // Adds an edge to an undirected graph
    void addEdge(int source, int destination)
    {
        adjacencylist[source].insert(destination);
        adjacencylist[destination].insert(source);
    }

    void printGraph()
    {
        for (int i = 0; i < vertices; i++)
        {
            unordered_set<int> neighbours = adjacencylist[i];
            cout << endl
                 << i << " : ";

            unordered_set<int>::iterator itr;
            for (itr = neighbours.begin(); itr != neighbours.end(); itr++)
            {
                cout << *itr << " ";
            }
            cout << endl;
        }
    }

    // Searches for a given edge in the graph
    bool isEdge(int source, int destination)
    {
        unordered_set<int>::iterator itr;
        itr = adjacencylist[source].find(destination);
        if (itr == adjacencylist[source].end())
            return false;
        else
            return true;
    }

    int getDegree(int vertex)
    {
        return adjacencylist[vertex].size();
    }

    int DSaturGraphColoring()
    {
        int maximum_degree_node = -1;
        int max_degree = 0;

        set<int> allColorsSET;

        for (int i = 0; i < vertices; i++)
        {
            timeslots[i] = -1;
            if (getDegree(i) > max_degree)
            {
                max_degree = getDegree(i);
                maximum_degree_node = i;
            }
        }
        //cout << "max dgree "<<maximum_degree_course << " : " << max_degree << endl;

        //Select a vertex of maximal degree and colour it with the first colour.
        int color = 0;
        timeslots[maximum_degree_node] = color;
        //dummy node isolated
        timeslots[0] = color;
        allColorsSET.insert(color);

        int nodeToBeColoured = vertices - 1;

        while (nodeToBeColoured > 0)
        {
            //find a vertex with the highest degree of saturation. Further ties are broken randomly.
            int max_saturation_degree = 0;
            int max_saturation_node = -1;
            for (int i = 1; i < vertices; i++)
            {
                if (timeslots[i] == -1)
                {
                    set<int> tempNeighboursColorsSET;
                    unordered_set<int>::iterator itr;
                    for (itr = adjacencylist[i].begin(); itr != adjacencylist[i].end(); itr++)
                    {
                        if (timeslots[*itr] != -1)
                        {
                            tempNeighboursColorsSET.insert(timeslots[*itr]);
                        }
                    }

                    if (max_saturation_degree == tempNeighboursColorsSET.size())
                    {
                        //Break ties by considering that vertex with the highest degree.
                        if (getDegree(i) > getDegree(max_saturation_node))
                        {
                            max_saturation_degree = tempNeighboursColorsSET.size();
                            max_saturation_node = i;
                        }
                        else if (getDegree(i) == getDegree(max_saturation_node))
                        {
                            int randomOfTwo = getRandomBetweenTwo(i,max_saturation_node);
                            max_saturation_degree = tempNeighboursColorsSET.size();
                            max_saturation_node = randomOfTwo;
                        }
                    }
                    else if (max_saturation_degree < tempNeighboursColorsSET.size())
                    {
                        max_saturation_degree = tempNeighboursColorsSET.size();
                        max_saturation_node = i;
                    }
                }
            }

            if (max_saturation_node == -1)
            {
                break;
            }

            //cout<<"yo5"<<endl;

            //Loop through the colour classes created so far, and colour the selected vertex with the first suitable colour.
            set<int> availableColorsSet = allColorsSET;

            set<int> neighboursColorsSET;
            unordered_set<int>::iterator itr;

            //cout<<"yo6"<<endl;
            for (itr = adjacencylist[max_saturation_node].begin(); itr != adjacencylist[max_saturation_node].end(); itr++)
            {
                //cout<<"sop : "<<max_saturation_node<<endl;
                if (timeslots[*itr] != -1)
                {
                    //cout<<"sop2"<<endl;
                    neighboursColorsSET.insert(timeslots[*itr]);
                }
            }

            set<int>::iterator sitr;
            for (sitr = neighboursColorsSET.begin(); sitr != neighboursColorsSET.end(); sitr++)
            {
                availableColorsSet.erase(*sitr);
            }

            if (availableColorsSet.empty())
            {
                color += 1;
                timeslots[max_saturation_node] = color;
                nodeToBeColoured--;
                allColorsSET.insert(color);
            }
            else
            {
                sitr = availableColorsSet.begin();
                timeslots[max_saturation_node] = *sitr;
                nodeToBeColoured--;
            }
        }

        /*for (int i = 0; i < vertices; i++)
        {
            cout << "Course : " << i << " ----->  Timeslot : " << timeslots[i] << endl;
        }*/
        //cout<<"Timeslots : "<<allColorsSET.size()<<endl;
        return allColorsSET.size();

        //tuple <int,int> tp;
    }

    int largestDegreeGraphColoring()
    {
        //timeslots[0] = 0;

        //pair(degree,course)
        priority_queue<pair<int, int>> pq;
        pair<int, int> top;

        for (int i = 0; i < vertices; i++)
        {
            pq.push(make_pair(getDegree(i), i));
        }

        //this array tracks if any of the node has same color in the adjacency vertices
        bool availableAdjacent[vertices];
        for (int i = 0; i < vertices; i++)
        {
            availableAdjacent[i] = false;
            timeslots[i] = -1;
        }
        //availableAdjacent[0] = true;

        while (!pq.empty())
        {
            top = pq.top();
            //cout << top.first << " " << top.second<<endl;
            int course_name = top.second;

            unordered_set<int>::iterator itr;

            for (itr = adjacencylist[course_name].begin(); itr != adjacencylist[course_name].end(); itr++)
            {
                if (timeslots[*itr] != -1)
                {
                    availableAdjacent[timeslots[*itr]] = true;
                }
            }

            int color;
            for (color = 0; color < vertices; color++)
            {
                if (availableAdjacent[color] == false)
                {
                    break;
                }
            }

            timeslots[course_name] = color;

            for (itr = adjacencylist[course_name].begin(); itr != adjacencylist[course_name].end(); itr++)
            {
                if (timeslots[*itr] != -1)
                {
                    availableAdjacent[timeslots[*itr]] = false;
                }
            }

            pq.pop();
        }

        int total_timeslot = 0;
        for (int i = 0; i < vertices; i++)
        {
            if (total_timeslot < timeslots[i])
            {
                total_timeslot = timeslots[i];
            }
            //cout << "Course : " << i << " ----->  Timeslot : " << timeslots[i] << endl;
        }

        //cout<<"Timeslots : "<<total_timeslot<<endl;
        return total_timeslot + 1;
    }

    void kempeChain(int r1, int r2)
    {
        vector<int> kempeChainNodes;
        queue<int> q;

        //cout<<r1<<"  "<<r2<<endl;

        for (int i = 0; i < vertices; i++)
        {
            discovered[i] = false;
        }

        int color1 = timeslots[r1];
        int color2 = timeslots[r2];

        discovered[r1] = true;
        discovered[r2] = true;

        kempeChainNodes.push_back(r1);
        kempeChainNodes.push_back(r2);

        q.push(r1);
        q.push(r2);

        while (!q.empty())
        {
            int f = q.front();
            q.pop();

            unordered_set<int>::iterator uitr;
            for (uitr = adjacencylist[f].begin(); uitr != adjacencylist[f].end(); uitr++)
            {
                if (discovered[*uitr] == false && (timeslots[*uitr] == color1 || timeslots[*uitr] == color2))
                {
                    q.push(*uitr);
                    discovered[*uitr] = true;
                    kempeChainNodes.push_back(*uitr);
                }
            }
        }

        kempeNodes = kempeChainNodes;
    }

    double penanltyCalculation(vector<vector<string>> &vec)
    {
        //w0=16, w1=8, w2=4, w3=2 and w4=1
        int penalty = 0;
        for (int i = 0; i < vec.size(); i++)
        {
            for (int j = 0; j < vec[i].size() - 1; j++)
            {
                for (int k = j + 1; k < vec[i].size(); k++)
                {
                    penalty += pow(2, 5 - abs(timeslots[stoi(vec[i][j])] - timeslots[stoi(vec[i][k])]));
                }
            }
        }

        double average_penalty = (double)penalty / vec.size();

        //cout<<"Penalty : "<<average_penalty<<endl;

        return average_penalty;
    }

    ~Graph()
    {
    }
};

// Driver code
int main()
{
    //extract course numbers from .crs file
    string s;
    int linecount;

    vector<pair<string, string>> filenames;
    filenames.push_back(make_pair("car-s-91.crs", "car-s-91.stu"));
    filenames.push_back(make_pair("car-f-92.crs", "car-f-92.stu"));
    filenames.push_back(make_pair("kfu-s-93.crs", "kfu-s-93.stu"));
    filenames.push_back(make_pair("tre-s-92.crs", "tre-s-92.stu"));
    filenames.push_back(make_pair("yor-f-83.crs", "yor-f-83.stu"));

    vector<string> solutionFiles;
    solutionFiles.push_back("car91.sol");
    solutionFiles.push_back("car92.sol");
    solutionFiles.push_back("kfu93.sol");
    solutionFiles.push_back("tre92.sol");
    solutionFiles.push_back("yor83.sol");

    cout << "                                  Largest Degree                                                         Desatur" << endl;
    cout << "Filename       Timeslots      Penalty(before kempe)       Penalty(after kempe)       Timeslots      Penalty(before kempe)       Penalty(after kempe)" << endl;

    for (int i = 0; i < 5; i++)
    {
        ifstream crs_in(filenames[i].first);
        while (getline(crs_in, s))
        {
            if (!s.empty())
                linecount++;
        }

        crs_in.close();

        int numOfCourses = linecount;

        Graph *g = new Graph(numOfCourses + 1);

        //graph creation from .stu file
        //line by line read

        ifstream infile(filenames[i].second);

        string line;

        //for storing file line
        vector<vector<string>> each_student_courses_vector;

        while (getline(infile, line))
        {
            if (!line.empty())
            {
                //vector<string> each_line;
                istringstream iss(line);

                vector<string> each_line{istream_iterator<string>{iss}, istream_iterator<string>{}};

                each_student_courses_vector.push_back(each_line);

                int len = each_line.size();

                for (int i = 0; i < len - 1; i++)
                {
                    for (int j = i + 1; j < len; j++)
                    {
                        g->addEdge(stoi(each_line[i]), stoi(each_line[j]));
                    }
                }
            }
        }

        infile.close();

        //g->printGraph();

        //real execution

        /******************************************/
        //cout<<"Largest Degree "<<endl;
        int colors_LD = g->largestDegreeGraphColoring();
        //cout<<"Timeslots : "<<colors_LD<<endl;
        double avg_penalty_LD = g->penanltyCalculation(each_student_courses_vector);
        double before_LD = avg_penalty_LD;
        //cout << "Average penalty(without kempe) : " << avg_penalty_LD << endl;

        for (int i = 0; i < 5000; i++)
        {
            int random_node = *select_randomly(g->allVerticesVector.begin(), g->allVerticesVector.end());
            int color1 = g->timeslots[random_node];

            vector<int> random_node_nighbours;
            unordered_set<int>::iterator uitr;
            for (uitr = g->adjacencylist[random_node].begin(); uitr != g->adjacencylist[random_node].end(); uitr++)
            {
                random_node_nighbours.push_back(*uitr);
            }

            if (random_node_nighbours.size() == 0)
                continue;

            //pick an adjacent neighbour
            int random_neighbour = *select_randomly(random_node_nighbours.begin(), random_node_nighbours.end());
            int color2 = g->timeslots[random_neighbour];

            kempeNodes.clear();

            g->kempeChain(random_node, random_neighbour);

            vector<int>::iterator vitr;
            //color swap
            for (vitr = kempeNodes.begin(); vitr != kempeNodes.end(); vitr++)
            {
                if (g->timeslots[*vitr] == color1)
                {
                    g->timeslots[*vitr] = color2;
                }
                else if (g->timeslots[*vitr] == color2)
                {
                    g->timeslots[*vitr] = color1;
                }
            }

            double temp_penalty = g->penanltyCalculation(each_student_courses_vector);

            //revert colors if penalty is larger
            if ((double)temp_penalty >= avg_penalty_LD)
            {
                for (vitr = kempeNodes.begin(); vitr != kempeNodes.end(); vitr++)
                {
                    if (g->timeslots[*vitr] == color1)
                    {
                        g->timeslots[*vitr] = color2;
                    }
                    else if (g->timeslots[*vitr] == color2)
                    {
                        g->timeslots[*vitr] = color1;
                    }
                }
            }
            else
            {
                avg_penalty_LD = temp_penalty;
                //cout << "Average penalty : " << avg_penalty << endl;
            }
        }
        //cout << "Average penalty( with kempe ) : " << avg_penalty_LD << endl;

        /************************************/

        //cout<<"Desatur "<<endl;
        int colors_d = g->DSaturGraphColoring();
        double avg_penalty_Dsatur = g->penanltyCalculation(each_student_courses_vector);
        double before_d = avg_penalty_Dsatur;
        //cout << "Average penalty : " << avg_penalty_Dsatur << endl;

        //loop for 1000 times
        for (int i = 0; i < 5000; i++)
        {
            int random_node = *select_randomly(g->allVerticesVector.begin(), g->allVerticesVector.end());
            int color1 = g->timeslots[random_node];

            vector<int> random_node_nighbours;
            unordered_set<int>::iterator uitr;
            for (uitr = g->adjacencylist[random_node].begin(); uitr != g->adjacencylist[random_node].end(); uitr++)
            {
                random_node_nighbours.push_back(*uitr);
            }

            if (random_node_nighbours.size() == 0)
                continue;

            //pick an adjacent neighbour
            int random_neighbour = *select_randomly(random_node_nighbours.begin(), random_node_nighbours.end());
            int color2 = g->timeslots[random_neighbour];

            kempeNodes.clear();

            g->kempeChain(random_node, random_neighbour);

            vector<int>::iterator vitr;
            //color swap
            for (vitr = kempeNodes.begin(); vitr != kempeNodes.end(); vitr++)
            {
                if (g->timeslots[*vitr] == color1)
                {
                    g->timeslots[*vitr] = color2;
                }
                else if (g->timeslots[*vitr] == color2)
                {
                    g->timeslots[*vitr] = color1;
                }
            }

            double temp_penalty = g->penanltyCalculation(each_student_courses_vector);

            //revert colors if penalty is larger
            if ((double)temp_penalty >= avg_penalty_Dsatur)
            {
                for (vitr = kempeNodes.begin(); vitr != kempeNodes.end(); vitr++)
                {
                    if (g->timeslots[*vitr] == color1)
                    {
                        g->timeslots[*vitr] = color2;
                    }
                    else if (g->timeslots[*vitr] == color2)
                    {
                        g->timeslots[*vitr] = color1;
                    }
                }
            }
            else
            {
                avg_penalty_Dsatur = temp_penalty;
                //cout << "Average penalty : " << avg_penalty << endl;
            }
        }

        cout << filenames[i].second << "      " << colors_LD << "               " << before_LD << "                     " << avg_penalty_LD << "                   " << colors_d << "                  " << before_d << "                    " << avg_penalty_Dsatur << endl;


        //optional 
        
        //ouput in .soln file
        /*ofstream ofs(solutionFiles[i]);
        for (int i = 1; i < numOfCourses + 1; i++)
        {
            ofs << i << "  " << g->timeslots[i] << endl;
        }*/
    }
    return 0;
}
