//
// Not fully done yet!! 
// Need to: 
// add separate lives, scores, change color of boxes, not end game when one person dies.. 
#include <iostream> 
using namespace std; 
#include <allegro5\allegro.h>
#include <allegro5\allegro_primitives.h>
#include "jump3.h"
#include <allegro5\allegro_font.h> 
#include <allegro5\allegro_ttf.h> 

//GLOBALS ================================
const int GRAVITY = -1;
const int MAX_SPEED = 15; // Keeping speed positive, but gravity negative because of physics
const int WIDTH = 400;
const int HEIGHT = 600;
const int NUM_LEDGES = 100; // maximum number of ledges on the screen
const int NUM_BULLETS = 15; 
const int JUMP_RANGE = 50;	// Max displacement between 2 ledges. 
// however, currently, the max displacement of the block is just an estimate that's >= 100. 
// keys are for 2 players. 
// (up, left, right, shoot left, shoot right)
enum KEYS{UP, LEFT, RIGHT,K, L, W, A, D, C, V};
bool keys[10] = {false, false, false, false, false, false, false, false, false, false};
bool landing = false; // if landing is true, it means the block needs to land. 
int score = 0; 
const int SCREEN_SPEED = 1; 

//prototypes =============================
void InitBlock(Block &block);
void DrawBlock(Block &block);
void MoveBlockUp(Block &block);
void MoveBlockLeft(Block &block);
void MoveBlockRight(Block &block);
void UpdateBlock(Block &block, Ledge ledges[]);	// Updates speed in y direction
inline void StopBlock(Block &block)
{
	block.speedy = 0; 
}

void InitLedge(Ledge ledges[], int size);
void DrawLedge(Ledge ledges[], int size);
void UpdateLedge(Ledge ledges[], int size);
bool LandLedge(Ledge ledges[], int lsize, Block &block);
int prev_ledge_height = HEIGHT;

void InitBullet( Bullet bullets[], int size); 
void DrawBullet( Bullet bullets[], int size); // draw the bullet to screen
void UpdateBullet( Bullet bullets[], int size); // update bullet to screen
void FireBulletR(Bullet bullets[], int size, Block &block);
void FireBulletL(Bullet bullets[], int size, Block &block);
void CollideBullet(Bullet bullet[], int bSize, Block &block);

void UpdateScroll(Block &block,Block &block2, Ledge ledges[], int Lsize);


int main(void)
{
	//primitive variables
	bool done = false;
	bool redraw = true;
	bool isGameOver = false;
	const int FPS = 60;

	//object variables
	Block block, block2; // initialize block
	Ledge ledges[NUM_LEDGES]; // initialize an array of ledges. 
	Bullet bullets[NUM_BULLETS]; 

	//Allegro variables
	ALLEGRO_DISPLAY *display = NULL;
	ALLEGRO_EVENT_QUEUE *event_queue = NULL;
	ALLEGRO_TIMER *timer = NULL;
	ALLEGRO_FONT *font18 = NULL; 
	//Initialization Functions
	if(!al_init())
		return -1;
	display = al_create_display(WIDTH, HEIGHT);
	if(!display)
		return -1;

	// install primitives and keyboards
	al_init_primitives_addon();
	al_install_keyboard();
	al_init_font_addon();
	al_init_ttf_addon(); 
	
	event_queue = al_create_event_queue();
	timer = al_create_timer(1.0/FPS);

	srand(time(NULL));
	// initialize block and ledges
	InitBlock(block);
	InitBlock(block2); 
	InitLedge(ledges, NUM_LEDGES);
	InitBullet(bullets, NUM_BULLETS); 
	font18 = al_load_font("Arial.ttf", 18, 0); 
	// register keyboard, timer and display. 
	al_register_event_source(event_queue, al_get_keyboard_event_source());
	al_register_event_source(event_queue, al_get_timer_event_source(timer));
	al_register_event_source(event_queue, al_get_display_event_source(display));

	al_start_timer(timer);
	while(!done)
	{
		ALLEGRO_EVENT ev;
		al_wait_for_event(event_queue, &ev);
		// when it's time. 
		if(ev.type == ALLEGRO_EVENT_TIMER)
		{
			redraw = true;
			if( !landing&& block2.y >= (HEIGHT - block2.boundy)) //Start upwards motion if block above floor
																		// during initialization. 
					block.y -=SCREEN_SPEED;
			if( !landing&& block2.y >= (HEIGHT - block2.boundy)) //Start upwards motion if block above floor
																		// during initialization.
					block2.y -=SCREEN_SPEED;
			if(!landing &&keys[UP])
			{
				MoveBlockUp(block);   
				MoveBlockUp(block2);
			}
			if(keys[LEFT])
				MoveBlockLeft(block);
			if(keys[RIGHT])
				MoveBlockRight(block);
			// if block is not dead after all the movement. 
			if(keys[K])
				FireBulletL(bullets, NUM_BULLETS, block); 
			if(keys[L])
				FireBulletR(bullets, NUM_BULLETS, block);
			if(keys[A])
				MoveBlockLeft(block2);
			if(keys[D])
				MoveBlockRight(block2);
			// if block is not dead after all the movement. 
			if(keys[C])
				FireBulletL(bullets, NUM_BULLETS, block2); 
			if(keys[V])
				FireBulletR(bullets, NUM_BULLETS, block2);
			if(!isGameOver)
			{
				if(LandLedge(ledges, NUM_LEDGES, block)) 
				{
					StopBlock(block);
					StopBlock(block2); 
					if (keys[UP])
					MoveBlockUp(block);
					if (keys[W])
					MoveBlockUp(block2); 
				}
				UpdateBlock(block, ledges);
				UpdateBlock(block2, ledges); 
				UpdateScroll(block, block2, ledges, NUM_LEDGES); 
				UpdateLedge(ledges, NUM_LEDGES); 
				UpdateBullet(bullets, NUM_BULLETS); 
				CollideBullet(bullets, NUM_BULLETS, block);
				CollideBullet(bullets, NUM_BULLETS, block2);
				if(block.y - 2 * block.boundy > HEIGHT)        //If top edge of block is below the screen
				{
					isGameOver = true;
					landing = false;  
				}
			}
		}

		// if user closes the window by mouse. 
		else if(ev.type == ALLEGRO_EVENT_DISPLAY_CLOSE)
			done = true;


		else if(ev.type == ALLEGRO_EVENT_KEY_DOWN)
		{
			switch(ev.keyboard.keycode)
			{	
				case ALLEGRO_KEY_ESCAPE:
					done = true;
					break;
				case ALLEGRO_KEY_UP:								
					keys[UP] = true;
					break;
				case ALLEGRO_KEY_LEFT:
					keys[LEFT] = true;
					break;
				case ALLEGRO_KEY_RIGHT:
					keys[RIGHT] = true;
					break;
				case ALLEGRO_KEY_K:
					keys[K] = true;
					break; 
				case ALLEGRO_KEY_L:
					keys[L] = true; 
					break; 
				case ALLEGRO_KEY_W:								
					keys[W] = true;
					break;
				case ALLEGRO_KEY_A:
					keys[A] = true;
					break;
				case ALLEGRO_KEY_D:
					keys[D] = true;
					break;
				case ALLEGRO_KEY_C:
					keys[C] = true;
					break; 
				case ALLEGRO_KEY_V:
					keys[V] = true; 
					break; 
			}
		}

		else if(ev.type == ALLEGRO_EVENT_KEY_UP)
		{
			switch(ev.keyboard.keycode)
			{
				case ALLEGRO_KEY_UP:								
					keys[UP] = false;
					break;
				case ALLEGRO_KEY_LEFT:
					keys[LEFT] = false;
					break;
				case ALLEGRO_KEY_RIGHT:
					keys[RIGHT] = false;
					break;
				case ALLEGRO_KEY_K:
					keys[K] = false;
					break; 
				case ALLEGRO_KEY_L:
					keys[L] = false; 
					break; 
				case ALLEGRO_KEY_W:								
					keys[W] = false;
					break;
				case ALLEGRO_KEY_A:
					keys[A] = false;
					break;
				case ALLEGRO_KEY_D:
					keys[D] = false;
					break;
				case ALLEGRO_KEY_C:
					keys[C] = false;
					break; 
				case ALLEGRO_KEY_V:
					keys[V] = false; 
					break; 
			}
		}
		//Check status of block each time through the while loop
		if(redraw && al_is_event_queue_empty(event_queue))
		{
			redraw = false;
			if(!isGameOver)
			{	
				DrawBlock(block);
				DrawBlock(block2); 
				DrawLedge(ledges, NUM_LEDGES);
				DrawBullet(bullets, NUM_BULLETS); 
			}
			else
			{
				// if game is over
			al_draw_textf(font18, al_map_rgb(255, 0, 0), WIDTH/2, HEIGHT/2, ALLEGRO_ALIGN_CENTRE, 
             "GAME OVER. SCORE: %i\n", score); 
			}
			al_flip_display();
			al_clear_to_color(al_map_rgb(0,0,0));
		}
	}
	al_destroy_display(display);
	return 0;
}

void InitBlock(Block &block)
{
	block.x = WIDTH / 2;
	block.y = HEIGHT - 30;	//Block edge will touch floor, not centre
	block.speedx = 7;
	block.speedy = 0;
	block.boundx = 10;
	block.boundy = 10;
	block.live = true;
}
void DrawBlock(Block &block)
{
	al_draw_filled_rectangle(block.x-block.boundx, block.y-block.boundy, 
		block.x + block.boundx, block.y + block.boundy, al_map_rgb(0,0,255));
}
void MoveBlockUp(Block &block)
{
	  landing = true; // start landing. 
	block.speedy = MAX_SPEED;
}
void MoveBlockLeft(Block &block)
{
	block.x -= block.speedx;
	if(block.x < 0)
		block.x = WIDTH; //Wraps around to other side of display
}
void MoveBlockRight(Block &block)
{
	block.x += block.speedx;								
	if(block.x > WIDTH)
		block.x = 0;										//Wraps around to other side of display
}
void UpdateBlock(Block &block, Ledge ledges[])	//DOES NOT MODIFY block.y AT ALL
{
	// if block is in the air, or just hitting floor
	if (landing) 
	{ //on it's way down.
		// if block.speedy > 0 : block is moving up 
		// if block.speedy = 0 : block is at it's highest point
		// if block.speedy < 0 : block is falling down
		block.y -= block.speedy;
		// block.speedy keeps adding gravity where gravity < 0. 
		block.speedy += GRAVITY;
	}
	if ( block.y <= 0)
		{
			block.y = 0;
			block.speedy = 0; 
		}
}

void InitLedge(Ledge ledges[], int size)
{
	// size is the number of ledges 
	for(int i = 0; i < size; i++)
	{
		ledges[i].x = rand() % WIDTH; // the number generated will be from 0 to width - 1
									// if you want number to include width, use ( rand() % (WIDTH +1) )
									// if you want to exclude 0, add 1. Use ( (rand() % (WIDTH)) + 1)  
		ledges[i].y = prev_ledge_height - ((rand() % (JUMP_RANGE/2)) + JUMP_RANGE/2); //ledges[i] will be within JUMP_RANGE above ledges[i-1]
		prev_ledge_height = ledges[i].y; // update previous ledge height. 
		ledges[i].boundx = 20;
		ledges[i].boundy = 20;
		ledges[i].live = true;
	}
}
void DrawLedge(Ledge ledges[], int size)
{
	for(int i = 0; i < size; i++)
	{
		// draw a rectangle for each ledge. 
		al_draw_filled_rectangle(ledges[i].x - ledges[i].boundx,
			ledges[i].y - ledges[i].boundy,
			ledges[i].x + ledges[i].boundx,
			ledges[i].y - 5,
			al_map_rgb(255, 255, 255));
	}
} 
bool LandLedge2(Ledge ledges[], int lsize, Block &block)
{
	int i = 0;
	// while block's height is higher than some ledges height. 
	while(block.y >= ledges[i++].y)
	{ // if it is falling
		if(block.speedy <= 0 && 
			// and its colliding with the current block vertically 
		((block.y + block.boundy) >= (ledges[i].y - ledges[i].boundy)) &&
		((block.y + block.boundy) <= (ledges[i].y + ledges[i].boundy)) &&	//Check if block is landing on a ledge
		// and colliding with current block horizontally. 
		(
			(((block.x - block.boundx <= (ledges[i].x + ledges[i].boundx)) && 
		 block.x-block.boundx >= (ledges[i].x - ledges[i].boundx)))    ||
			(((block.x + block.boundx <= (ledges[i].x + ledges[i].boundx)) &&
		 block.x+block.boundx >= (ledges[i].x - ledges[i].boundx)))
		 )
		)
		{
			return true;
		}
	}
	return false;
}

bool LandLedge(Ledge ledges[], int lsize, Block &block)
{
    for(int i = 0; i < lsize; i++)
    {
		// if block is falling. 
        if(block.speedy < 0)
			// if block is within the y bound
        {    if ((block.y + block.boundy >= ledges[i].y - ledges[i].boundy) &&
                (block.y + block.boundy <= ledges[i].y + ledges[i].boundy))
            { // if x is within the bound
                if ((block.x + block.boundx <= ledges[i].x + ledges[i].boundx) &&
                    (block.x + block.boundx >= ledges[i].x - ledges[i].boundx))
                {
                    cout<<"True!" <<endl;
                    return true;
                }
				// or if x is within the bound of the other one. 
                else if ((block.x - block.boundx >= ledges[i].x - ledges[i].boundx) &&
                    (block.x - block.boundx <= ledges[i].x + ledges[i].boundx))
                {
                    cout << "True!" << endl;
                return true;
                }
            }
        }
    }
    cout<<"FALSE"<<endl;
    return false;
} 

void UpdateLedge(Ledge ledges[], int size)
{
   // prev_ledge_height += HEIGHT/4; 
	for (int i = 0; i < size; i ++)
	{
		if ( ledges[i].live == false) 
		{
		ledges[i].x = rand() % WIDTH; // the number generated will be from 0 to width - 1
		// if you want number to include width, use ( rand() % (WIDTH +1) )
		// if you want to exclude 0, add 1. Use ( (rand() % (WIDTH)) + 1)  
		ledges[i].y = prev_ledge_height - (rand() % (JUMP_RANGE)); //ledges[i] will be within JUMP_RANGE above ledges[i-1]
		prev_ledge_height = ledges[i].y; // update previous ledge height. 
		ledges[i].boundx = 20;
		ledges[i].boundy = 20;
		ledges[i].live = true;
		}
	}
}

void UpdateScroll(Block &block,Block &block2, Ledge ledges[], int Lsize)
{
		for (int i = 0; i < Lsize; i++)
		{
			ledges[i].y += SCREEN_SPEED;
			if (ledges[i].y > HEIGHT)
			{
				ledges[i].live = false; 
			}
		}
		prev_ledge_height += SCREEN_SPEED; 
		block.y += SCREEN_SPEED;
		block2.y += SCREEN_SPEED; 
		if (landing) score += SCREEN_SPEED; 
}

void InitBullet( Bullet bullets[], int size)
{
	for ( int i = 0; i < size; i++)
	{
		bullets[i].speed = 10; 
		bullets[i].live = false; 
	}
}
void DrawBullet( Bullet bullets[], int size)
{
	// if the bullet is alive, draw it. 
	for (int i = 0; i < size; i++)
	{
		if(bullets[i].live) // draw the bullet at it's current position. 
			al_draw_filled_circle(bullets[i].x, bullets[i].y, 2, al_map_rgb(255,255, 255)); 
	}
}
void UpdateBullet( Bullet bullets[], int size)
{
	for (int i = 0; i < size; i++)
	{
		// if the bullet [i] is alive, 
		if(bullets[i].live)
		{
			//bullets[i].y += SCREEN_SPEED;
			if (bullets[i].ID == 1)
			{
				bullets[i].x += bullets[i].speed; // continue moving it at to the right.
				if(bullets[i].x > WIDTH) // until it reaches the end of the screen
					bullets[i].live = false;  // then it is no longer alive. 
			}
			else
			{
				bullets[i].x -= bullets[i].speed; // continue moving it to the left.  
				if(bullets[i].x < 0) // until it reaches the end of the screen
					bullets[i].live = false;  // then it is no longer alive. 
			}
		}
	}
}
void FireBulletR(Bullet bullets[], int size, Block &block)
{
	for (int i = 0; i < size; i++)
	{  // find 1 bullet that is not alive, 
		if(!bullets[i].live)
		{	// fire the bullet. 
			bullets[i].ID = 1; 
			bullets[i].x = block.x + 17; // the initial position is 17 away from the ship's position
			bullets[i].y = block.y; // and same height as ship. 
			bullets[i].live = true; // change that bullet to be alive
			break; // exit the for loop so it doesn't fire the next bullet. 
		}
	}
}
void FireBulletL(Bullet bullets[], int size, Block &block)
{
	for (int i = 0; i < size; i++)
	{  // find 1 bullet that is not alive, 
		if(!bullets[i].live)
		{	// fire the bullet. 
			bullets[i].ID = 2; 
			bullets[i].x = block.x - 17; // the initial position is 17 away from the ship's position
			bullets[i].y = block.y; // and same height as ship. 
			bullets[i].live = true; // change that bullet to be alive
			break; // exit the for loop so it doesn't fire the next bullet. 
		}
	}
}
void CollideBullet(Bullet bullet[], int bSize, Block &block)
{
	for ( int i = 0; i < bSize; i++)
	{
		if(bullet[i].live)
		{
			if(block.live && bullet[i].ID == 2)
			{		// if bullet is between x and boundx of ship2
					// and between y and boundy of ship2
					// there is a collision. 
				if( bullet[i].x > (block.x - block.boundx) && 
					bullet[i].x < (block.x + block.boundx) &&
					bullet[i].y > (block.y - block.boundy) && 
					bullet[i].y < (block.y + block.boundy)) 
				{ 
					bullet[i].live = false; 
					block.x-=3; // move block to the left. 
				}
			}
			if(block.live && bullet[i].ID == 1 )
			{
				if( bullet[i].x > (block.x - block.boundx) && 
					bullet[i].x < (block.x + block.boundx) &&
					bullet[i].y > (block.y - block.boundy) && 
					bullet[i].y < (block.y + block.boundy)) 
				{ 
					bullet[i].live = false; 
					block.x += 3;// move block to the right. 
				}
			}
		}
	}
}

// */
