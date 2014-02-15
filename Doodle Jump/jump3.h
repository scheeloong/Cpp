struct Block
{
	int x;
	int y;
	int speedx;
	int speedy;
	int boundx;
	int boundy;
	bool live; // to determine if block is alive, if block drops dead, it means block is dead
				// which means game over.. 
};

struct Ledge
{
	int x;
	int y;
	int boundx;
	int boundy;
	bool live; // to determine if current ledge is alive. 
};

struct Bullet
{
	int ID;
	int x;
	int y; 
	bool live; 
	int speed; 
};