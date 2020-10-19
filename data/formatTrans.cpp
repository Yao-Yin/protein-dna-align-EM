#include <iostream>
#include <string>

using namespace std;

int main() {
    string curr;
    string res;
    cin >> curr;
    cout << curr;
    while(cin >> curr) {
        res += curr;
    }
    for (int i = 0; i < res.size(); i ++) {
        if(i % 60 == 0) cout << endl;
        cout << res[i];
    }
    return 0;
}