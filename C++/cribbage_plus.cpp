///Cribbage Scorekeeper (Up to Four Players)
///6/14/20


#include <iostream>

using namespace std;

int main() {
  int max_hands  = 5000;
  int score1 = 0;
  int score2 = 0;
  int score3 = 0;
  int score4 = 0;
  int hscore1 = 0;
  int hscore2 = 0;
  int hscore3 = 0;
  int hscore4 = 0;
  int point_lay;
  int who;
  int extra_points;
  int kitty;
  string p1name;
  string p2name;
  string p3name;
  string p4name;
  int num_players;
  string  corresp_string;

  cout << "How many players?" << endl;
  cin >> num_players;
  cout << endl;
  if(num_players<3){
    cout << "Enter Player 1's Name" << endl;
    cin >> p1name;
    cout << endl;
    cout << "Enter Player 2's Name" << endl;
    cin >> p2name;
    cout << endl;
    corresp_string = "(" +p1name+"=1/"+p2name+"=2)";
  }
  if(num_players==3){
    cout << "Enter Player 1's Name" << endl;
    cin >> p1name;
    cout << endl;
    cout << "Enter Player 2's Name" << endl;
    cin >> p2name;
    cout << endl;
    cout << "Enter Player 3's Name" << endl;
    cin >> p3name;
    cout << endl;
    corresp_string = "(" +p1name+"=1/"+p2name+"=2/"+p3name+"=3)";
  }
  if(num_players == 4){
    cout << "Enter Player 1's Name" << endl;
    cin >> p1name;
    cout << endl;
    cout << "Enter Player 2's Name" << endl;
    cin >> p2name;
    cout << endl;
    cout << "Enter Player 3's Name" << endl;
    cin >> p3name;
    cout << endl;
    cout << "Enter Player 4's Name" << endl;
    cin >> p4name;
    cout << endl;
    corresp_string = "(" +p1name+"=1/"+p2name+"=2/"+p3name+"=3/"+p4name+"=4)";
  }

  cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" <<endl;


    for(int ii = 0; ii<max_hands; ii++){
      for(int jj =0; jj<max_hands; jj++){
        cout << "Was there a point during play? (1=yes/0=no)" << endl;
        cin >> point_lay;
        if(point_lay == 1){
          cout << "Who got it? " <<corresp_string<< endl;
          cin >> who;
          if(who == 1){
            cout << "Enter points" << endl;
            cin>>extra_points;
            score1+=extra_points;
            cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" <<endl;
            if(score1>120){
              cout <<p1name <<" wins!" << endl;
              exit(1);
            }
          }else if(who == 2){
            cout << "Enter points" << endl;
            cin>>extra_points;
            score2+=extra_points;
            cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" <<endl;
            if(score2>120){
              cout << p2name <<" wins!" << endl;
              exit(1);
            }
          }else if(who == 3){
            if(who>num_players){
              cout << "There aren't that many players! Please try again.";
            }else{
              cout << "Enter points" << endl;
              cin>>extra_points;
              score3+=extra_points;
              cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" <<endl;
              if(score3>120){
                cout << p3name <<" wins!" << endl;
                exit(1);
              }
            }
          }else{
            if(who>num_players){
              cout << "There aren't that many players! Please try again.";
            }else{
              cout << "Enter points" << endl;
              cin>>extra_points;
              score4+=extra_points;
              cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" <<endl;
              if(score4>120){
                cout << p4name <<" wins!" << endl;
                exit(1);
              }
            }
          }
        }else{
          break;
        }
      }
      if(ii % num_players == 0){
        cout <<"What was "<< p1name <<"'s score this hand?" << endl;
        cin >> hscore1;
        score1+=hscore1;
        cout << p1name<<"'s score is " << score1 << endl;
        cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" <<endl;
        if(score1>120){
          cout << p1name << " wins!" << endl;
          exit(1);
        }
        cout <<"What was "<< p2name <<"'s score this hand?" << endl;
        cin >> hscore2;
        score2+=hscore2;
        if(num_players<3){
          cout << "What's the kitty?" << endl;
          cin >> kitty;
          score2+=kitty;
        }
        cout << p2name<<"'s score is " << score2 << endl;
        cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" <<endl;
        if(score2>120){
          cout << p2name << " wins!" << endl;
          exit(1);
        }
        if(num_players>2){
          cout <<"What was "<< p3name <<"'s score this hand?" << endl;
          cin >> hscore3;
          score3+=hscore3;
          if(num_players == 3){
            cout << "What's the kitty?" << endl;
            cin >> kitty;
            score3+=kitty;
          }
          cout << p3name<<"'s score is " << score3 << endl;
          cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" <<endl;
          if(score3>120){
            cout << p3name << " wins!" << endl;
            exit(1);
          }
        }
        if(num_players == 4){
          cout <<"What was "<< p4name <<"'s score this hand?" << endl;
          cin >> hscore4;
          score4+=hscore4;
          cout << "What's the kitty?" << endl;
          cin >> kitty;
          score4+=kitty;
          cout << p4name<<"'s score is " << score4 << endl;
          cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" <<endl;
          if(score4>120){
            cout << p4name << " wins!" << endl;
            exit(1);
          }
        }

      }
      if(ii % num_players == 1){
        cout <<"What was "<< p2name <<"'s score this hand?" << endl;
        cin >> hscore2;
        score2+=hscore2;
        cout << p2name<<"'s score is " << score2 << endl;
        cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" <<endl;
        if(score2>120){
          cout << p2name << " wins!" << endl;
          exit(1);
        }
        if(num_players>2){
          cout <<"What was "<< p3name <<"'s score this hand?" << endl;
          cin >> hscore3;
          score3+=hscore3;
          cout << p3name<<"'s score is " << score3 << endl;
          cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" <<endl;
          if(score3>120){
            cout << p3name << " wins!" << endl;
            exit(1);
          }
        }
        if(num_players == 4){
          cout <<"What was "<< p4name <<"'s score this hand?" << endl;
          cin >> hscore4;
          score4+=hscore4;
          cout << p4name<<"'s score is " << score4 << endl;
          cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" <<endl;
          if(score4>120){
            cout << p4name << " wins!" << endl;
            exit(1);
          }
        }
        cout <<"What was "<< p1name <<"'s score this hand?" << endl;
        cin >> hscore1;
        score1+=hscore1;
        cout << "What's the kitty?" << endl;
        cin >> kitty;
        score1+=kitty;
        cout << p1name<<"'s score is " << score1 << endl;
        cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" <<endl;
        if(score1>120){
          cout << p1name << " wins!" << endl;
          exit(1);
        }
      }
      if(num_players >2){
        if(ii % num_players == 2){
          if(num_players>2){
            cout <<"What was "<< p3name <<"'s score this hand?" << endl;
            cin >> hscore3;
            score3+=hscore3;
            cout << p3name<<"'s score is " << score3 << endl;
            cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" <<endl;
            if(score3>120){
              cout << p3name << " wins!" << endl;
              exit(1);
            }
          }
          if(num_players == 4){
            cout <<"What was "<< p4name <<"'s score this hand?" << endl;
            cin >> hscore4;
            score4+=hscore4;
            cout << p4name<<"'s score is " << score4 << endl;
            cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" <<endl;
            if(score4>120){
              cout << p4name << " wins!" << endl;
              exit(1);
            }
          }
          cout <<"What was "<< p1name <<"'s score this hand?" << endl;
          cin >> hscore1;
          score1+=hscore1;
          cout << p1name<<"'s score is " << score1 << endl;
          cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" <<endl;
          if(score1>120){
            cout << p1name << " wins!" << endl;
            exit(1);
          }
          cout <<"What was "<< p2name <<"'s score this hand?" << endl;
          cin >> hscore2;
          score2+=hscore2;
          cout << "What's the kitty?" << endl;
          cin >> kitty;
          score2+=kitty;
          cout << p2name<<"'s score is " << score2 << endl;
          cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" <<endl;
          if(score2>120){
            cout << p2name << " wins!" << endl;
            exit(1);
          }
        }
      }
      if(num_players == 4){
        if(ii % num_players == 3){
          if(num_players == 4){
            cout <<"What was "<< p4name <<"'s score this hand?" << endl;
            cin >> hscore4;
            score4+=hscore4;
            cout << p4name<<"'s score is " << score4 << endl;
            cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" <<endl;
            if(score4>120){
              cout << p4name << " wins!" << endl;
              exit(1);
            }
          }
          cout <<"What was "<< p1name <<"'s score this hand?" << endl;
          cin >> hscore1;
          score1+=hscore1;
          cout << p1name<<"'s score is " << score1 << endl;
          cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" <<endl;
          if(score1>120){
            cout << p1name << " wins!" << endl;
            exit(1);
          }
          cout <<"What was "<< p2name <<"'s score this hand?" << endl;
          cin >> hscore2;
          score2+=hscore2;
          cout << p2name<<"'s score is " << score2 << endl;
          cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" <<endl;
          if(score2>120){
            cout << p2name << " wins!" << endl;
            exit(1);
          }
          if(num_players>2){
            cout <<"What was "<< p3name <<"'s score this hand?" << endl;
            cin >> hscore3;
            score3+=hscore3;
            cout << "What's the kitty?" << endl;
            cin >> kitty;
            score3+=kitty;
            cout << p3name<<"'s score is " << score3 << endl;
            cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" <<endl;
            if(score3>120){
              cout << p3name << " wins!" << endl;
              exit(1);
            }
          }
        }
      }
    }
  return 0;
}
