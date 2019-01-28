// Object IDS
// for every object to have a unique ID
enum IDS{ PLAYER, BULLET, PLAYER2}; 
//          0 ,    1,       2
// Our Player
struct SpaceShip
{
	int ID; // its ID
	int x;  // position x
	int y;  // position y
	int lives; // current live
	int speed; // speed
	int boundx;  // bounds are used for collision. 
	int boundy;  // it determines if two objects collide. 
	int score;  // score
};

struct Bullet
{
	int ID; // its ID
	int x; // position x
	int y; // position y
	bool live; // current live
	int speed; // speed
};


struct Comet
{
	int ID;
	int x;
	int y; 
	bool live; 
	int speed;
	// bounds for collision detection
	int boundx;
	int boundy; 
}; 