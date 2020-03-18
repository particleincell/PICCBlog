#include <LiquidCrystal.h>

// specify pins to use
const int pinButton = 22;  
const int pinRed = 53;  
const int pinGreen = 51;  
const int pinBlue = 49;
const int rs = 12, en = 11, d4 = 5, d5 = 4, d6 = 3, d7 = 2;
LiquidCrystal lcd(rs, en, d4, d5, d6, d7);  // lcd object

void setup() {
  // put your setup code here, to run once:
  Serial.begin(57600);
  pinMode(pinRed,OUTPUT);
  pinMode(pinGreen,OUTPUT);
  pinMode(pinBlue,OUTPUT);
  lcd.begin(16, 2);    // set up LCD dimensions
  setLedColor(0,0,0);  // turn off the LED

  reset();    // initialize counters
}

int button_old = 0;    // the last state of the button
long N_in = 0;
long N_tot = 0;

void loop() {
  // put your main code here, to run repeatedly:
  int button = digitalRead(pinButton);
  if (!button && button_old) reset();
  button_old = button;

  // pick random point
  float x = random(RAND_MAX)/(float)RAND_MAX;
  float y = random(RAND_MAX)/(float)RAND_MAX;
  N_tot++;

  if (x*x+y*y<=1) {      // if inside
    setLedColor(0,1,0);  // green
    N_in++;
  }
  else {
    setLedColor(1,0,0);  // red
  }

  // update pi estimate on the LCD
  lcd.setCursor(9,0);
  lcd.print(N_tot);
  lcd.setCursor(4,1);
  float pi = 4*N_in/(float)N_tot;
  lcd.print(pi,4);
   
  delay(50);  // wait 50 ms to slow down the LED blinking rate
}

// illuminates the LED to a particular combination of r/g/b
void setLedColor(int r,int g, int b) {
 digitalWrite(pinRed,1-r);
 digitalWrite(pinGreen,1-g);
 digitalWrite(pinBlue,1-b);
}

void reset() {
 lcd.clear();
 lcd.setCursor(0,0);
 lcd.print("Samples:");
 lcd.setCursor(0,1);
 lcd.print("pi:");
 Serial.println("Resetting!");
 setLedColor(0,0,0);
 delay(200);
 setLedColor(0,0,1);
 delay(500);
 setLedColor(0,0,0);
 delay(200); 
 N_in = 0;
 N_tot = 0;
}
