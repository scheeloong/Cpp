// 
// First Game. SpaceShip (SideShooter)
#include <allegro5\allegro.h> 
#include <allegro5\allegro_primitives.h> 
#include <allegro5\allegro_font.h> // note program will crash if the specified font
#include <allegro5\allegro_ttf.h> // is not available. 
#include "Game2.h" 
#include <allegro5\allegro_image.h> 

// Globals
const int WIDTH = 800;
const int HEIGHT = 400;
const int NUM_BULLETS = 5; // maximum number of bullets.
// space for fire
enum KEYS{UP, DOWN, LEFT, RIGHT, SPACE, W, A, S, D,E};
bool keys[10] = {false, false, false, false, false, false, false, false, false, false}; // false == 0, true == 1; 

// function prototypes. 
void InitShip(SpaceShip &ship);  // initialize ship
void InitShip2(SpaceShip &ship2); // initialize ship2

void DrawShip(SpaceShip &ship, ALLEGRO_BITMAP * &image);  // draw ship to screen

void MoveShipUp(SpaceShip &ship);  // move the ship 
void MoveShipDown(SpaceShip &ship); 

void MoveShipLeft(SpaceShip &ship); 
void MoveShip2Left(SpaceShip &ship2);
void MoveShipRight(SpaceShip &ship);
void MoveShip2Right(SpaceShip &ship2);

void InitBullet(Bullet bullet[], int size); // initalize bullet with size
void DrawBullet(Bullet bullet[], int size); // draw bullet with size
void FireBullet(Bullet bullet[], int size, SpaceShip &ship); // fire the bullet from ship. 
void FireBullet2(Bullet bullet[], int size, SpaceShip &ship2);
void UpdateBullet(Bullet bullet[], int size); // updates the position of bullet, and if it still exists.
void UpdateBullet2(Bullet bullet[], int size);
// continuously updates the position of bullet until it reaches the end of screen . 
void CollideBullet(Bullet bullet[], int bSize, SpaceShip &ship2, SpaceShip &ship);
void CollideBullet2(Bullet bullet[], int b2Size, SpaceShip &ship, SpaceShip &ship2);



int main(void)
{
	// primitive variables
	bool done = false; // for game loop
	bool redraw = true;  // for redraw. 
	const int FPS = 60;  // Frames Per Second
	bool isGameOver = false; 
	// initialize
	// object variables
	SpaceShip ship, ship2;
	Bullet bullet[NUM_BULLETS]; 
	Bullet bullet2[NUM_BULLETS]; 
	// Initialize 
	ALLEGRO_DISPLAY *display = NULL; 
	ALLEGRO_EVENT_QUEUE *event_queue = NULL ;
	ALLEGRO_TIMER *timer = NULL;
	ALLEGRO_FONT *font18 = NULL; 
	ALLEGRO_BITMAP *image1 = NULL , *image2 = NULL; 
	if(!al_init()) return -1; 
	display = al_create_display(WIDTH, HEIGHT); 
		if(!display) return -1; 

	// initializing primitives (for drawing ship) 
	al_init_primitives_addon(); 
	al_install_keyboard(); 
	al_init_font_addon();
	al_init_ttf_addon(); 
	al_init_image_addon(); 
	event_queue = al_create_event_queue(); 
	timer = al_create_timer(1.0/FPS); 
	// initialize objects
	srand(time(NULL)); // sees number generator with current time. 
						// as time changes, the random number changes too. 
	InitShip(ship); 
	InitShip2(ship2); 
	InitBullet(bullet, NUM_BULLETS); // array of bullets for ship
	InitBullet(bullet2, NUM_BULLETS); // array of bullets for ship2
	font18 = al_load_font("Arial.ttf", 18, 0); 
	image1 = al_load_bitmap("Ship1.bmp"); 
	image2 = al_load_bitmap("Ship2.bmp"); 
	// register keyboard. 
	al_register_event_source(event_queue, al_get_keyboard_event_source()); 
	al_register_event_source(event_queue, al_get_timer_event_source(timer)); 
	al_register_event_source(event_queue, al_get_display_event_source(display)); 
	// START TIMER
	al_start_timer(timer); 
	while(!done)
	{
		ALLEGRO_EVENT ev; 
		// wait for event to come in from queue. 
		al_wait_for_event(event_queue, &ev); 
		// only happens once every 60 seconds max. 
		if (ev.type == ALLEGRO_EVENT_TIMER)
		{
			redraw = true; 
			// update the position of the ship, depending if the keys are pressed. 
			if(keys[UP]) MoveShipUp(ship2);
			if(keys[DOWN]) MoveShipDown(ship2);
			if(keys[LEFT]) MoveShip2Left(ship2);
			if(keys[RIGHT]) MoveShip2Right(ship2);
			if(keys[W]) MoveShipUp(ship);
			if(keys[S]) MoveShipDown(ship);
			if(keys[A]) MoveShipLeft(ship);
			if(keys[D]) MoveShipRight(ship);

			if(!isGameOver)
			{
				UpdateBullet(bullet, NUM_BULLETS);  // update the position of the bullet
				UpdateBullet2(bullet2, NUM_BULLETS); 
				// test if all the objects collide after all the objects move. 
				CollideBullet(bullet, NUM_BULLETS, ship2, ship); 
				CollideBullet2(bullet2, NUM_BULLETS, ship, ship2); 
				if(ship.lives <= 0 || ship2.lives <=0) isGameOver = true; 
			}
		}
		else if(ev.type == ALLEGRO_EVENT_DISPLAY_CLOSE) done = true; 
		else if (ev.type == ALLEGRO_EVENT_KEY_DOWN)
		{
			switch(ev.keyboard.keycode)
			{
			case ALLEGRO_KEY_ESCAPE:
				done = true;
				break;
			case ALLEGRO_KEY_UP:
				keys[UP] = true;
				break;
			case ALLEGRO_KEY_DOWN:
				keys[DOWN] = true;
				break;
			case ALLEGRO_KEY_LEFT:
				keys[LEFT] = true;
				break;
			case ALLEGRO_KEY_RIGHT:
				keys[RIGHT] = true;
				break;
			case ALLEGRO_KEY_SPACE:
				keys[SPACE] = true; 
				FireBullet2(bullet2, NUM_BULLETS, ship2); // break to fire one bullet only
				break; 
			case ALLEGRO_KEY_W:
				keys[W] = true;
				break;
			case ALLEGRO_KEY_A:
				keys[A] = true;
				break;
			case ALLEGRO_KEY_S:
				keys[S] = true;
				break;
			case ALLEGRO_KEY_D:
				keys[D] = true;
				break;
			case ALLEGRO_KEY_E:
				keys[E] = true; 
				FireBullet(bullet, NUM_BULLETS, ship); // break to fire one bullet only
				break; 
			}
		}

		else if(ev.type == ALLEGRO_EVENT_KEY_UP)
		{
			switch(ev.keyboard.keycode)
			{
				case ALLEGRO_KEY_ESCAPE:
					done = true;
					break;
				case ALLEGRO_KEY_UP:
					keys[UP] = false;
					break;
				case ALLEGRO_KEY_DOWN:
					keys[DOWN] = false;
					break;
				case ALLEGRO_KEY_LEFT:
					keys[LEFT] = false;
					break;
				case ALLEGRO_KEY_RIGHT:
					keys[RIGHT] = false;
					break;
				case ALLEGRO_KEY_SPACE:
					keys[SPACE] = false; 
					break; 
				case ALLEGRO_KEY_W:
					keys[W] = false;
					break;
				case ALLEGRO_KEY_A:
					keys[A] = false;
					break;
				case ALLEGRO_KEY_S:
					keys[S] = false;
					break;
				case ALLEGRO_KEY_D:
					keys[D] = false;
					break;
				case ALLEGRO_KEY_E:
					keys[E] = false; 
					// break to fire one bullet only
					break; 
			}
		}

		// if redraw is true and no queue for event. 
		if ( redraw && al_is_event_queue_empty(event_queue))
		{
			redraw = false; 

			if (!isGameOver)
			{
				// draws the latest position of bullet and ship that has been updated. 
				DrawShip(ship, image1); // Draw ship
				DrawShip(ship2, image2); 
				DrawBullet(bullet, NUM_BULLETS); // draw bullets 
				DrawBullet(bullet2, NUM_BULLETS); 
				al_draw_textf(font18, al_map_rgb(255,0,255), 5, 5, 0, "Player 1 Lives left: %i \nScore: %i\n", ship.lives, ship.score); 
				al_draw_textf(font18, al_map_rgb(255,0,255), 50, 50, 0, "Player 2Lives left: %i \nScore: %i\n", ship2.lives, ship2.score);
			}
			else 
			{ al_draw_textf(font18, al_map_rgb(0,255,255), WIDTH/2, HEIGHT/2, ALLEGRO_ALIGN_CENTRE,
							"GAME OVER. P1 Final Score: %i\n P2 Final Score %i",ship.score, ship2.score); 

			}
			al_flip_display(); 
			al_clear_to_color(al_map_rgb(0, 0, 0)); 
		}
	}
	al_destroy_display(display); 
	return 0;
}

// initializing ship function
void InitShip(SpaceShip &ship)
{
	ship.ID = PLAYER;
	ship.x = 20; // initial positions
	ship.y = HEIGHT / 2;
	ship.lives = 10; 
	ship.speed = 7; // speed that it moves. 
	ship.boundx = 18;
	ship.boundy = 18;
	ship.score = 0; 
}
// initializing ship2 function
void InitShip2(SpaceShip &ship2)
{
	ship2.ID = PLAYER2; 
	ship2.x = WIDTH - 20; 
	ship2.y = HEIGHT / 2; 
	ship2.lives = 10; 
	ship2.speed = 7; 
	ship2.boundx = 18; 
	ship2.boundy = 18; 
	ship2.score = 0;
}


void DrawShip(SpaceShip &ship, ALLEGRO_BITMAP * &image)
{
	if (ship.lives <= 0)
	{
		ship.y = -100;
	}

	// Drawing the ship. 
	else
	{	
		al_convert_mask_to_alpha(image, al_map_rgb(255,255,255)); 
		al_draw_bitmap(image, ship.x - al_get_bitmap_width(image)/2, ship.y - al_get_bitmap_height(image)/2, 0);
	}
}

void MoveShipUp(SpaceShip &ship)
{
	ship.y -= ship.speed; 
	if (ship.y < 0)  // if ship is at highest point
		ship.y = 0;  // it remains at highest point
}
void MoveShipDown(SpaceShip &ship)
{
	ship.y += ship.speed; 
	if (ship.y > HEIGHT) // if ship is at lowest point of window
		ship.y = HEIGHT;  // it remains at lowest point of window. 
}
void MoveShipLeft(SpaceShip &ship)
{
	ship.x -= ship.speed; 
	if (ship.x < 0) 
		ship.x = 0; // if ship is at furthest left of window, it remains there. 
}
void MoveShip2Left(SpaceShip &ship2)
{
	ship2.x -= ship2.speed; 
	if (ship2.x < WIDTH - 300) 
		ship2.x = WIDTH - 300; // if ship is at furthest left of window, it remains there. 
}
void MoveShipRight(SpaceShip &ship)
{
	ship.x += ship.speed; // move it by the amount of it's speed. 
	// limit how far player can go to the right. 
	if (ship.x > 300)  // if ship exceeds 300, 
		ship.x = 300; // it remains at 300. 
}
void MoveShip2Right(SpaceShip &ship2)
{
	ship2.x += ship2.speed; // move it by the amount of it's speed. 
	// limit how far player can go to the right. 
	if (ship2.x > WIDTH)  // if ship exceeds Width, 
		ship2.x = WIDTH; // it remains at Width.  
}

// initialize bullet for both ships. 
// movement of bullets is determined by FireBullet & FireBullet2
void InitBullet(Bullet bullet[], int size)
{
	// make number of bullets. 
	for (int i = 0; i < size; i++)
	{
		bullet[i].ID = BULLET; 
		bullet[i].speed = 10; 
		bullet[i].live = false; 
	}
}
// Draw the bullet
void DrawBullet(Bullet bullet[], int size)
{
	// if the bullet is alive, draw it. 
	for (int i = 0; i < size; i++)
	{
		if(bullet[i].live) // draw the bullet at it's current position. 
			al_draw_filled_circle(bullet[i].x, bullet[i].y, 2, al_map_rgb(255,255, 255)); 
	}
}

void FireBullet(Bullet bullet[], int size, SpaceShip &ship )
{
	for (int i = 0; i < size; i++)
	{  // find 1 bullet that is not alive, 
		if(!bullet[i].live)
		{	// fire the bullet. 
			bullet[i].x = ship.x + 17; // the initial position is 17 away from the ship's position
			bullet[i].y = ship.y; // and same height as ship. 
			bullet[i].live = true; // change that bullet to be alive
			break; // exit the for loop so it doesn't fire the next bullet. 
		}
	}
}
void FireBullet2(Bullet bullet[], int size, SpaceShip &ship2 )
{
	for (int i = 0; i < size; i ++)
	{  // find 1 bullet that is not alive, 
		if(!bullet[i].live)
		{	// fire the bullet. 
			bullet[i].x = ship2.x - 17; // the initial position is 17 away from the ship's position 
										// to the left. 
			bullet[i].y = ship2.y; // and same height as ship. 
			bullet[i].live = true; // change that bullet to be alive
			break; // exit the for loop so it doesn't fire the next bullet. 
		}
	}
}

// function that update position of bullet for ship1. 
void UpdateBullet(Bullet bullet[], int size)
{
	for (int i = 0; i < size; i++)
	{
		// if the bullet [i] is alive, 
		if(bullet[i].live)
		{
			bullet[i].x += bullet[i].speed; // continue moving it at to the right.  
			if(bullet[i].x > WIDTH) // until it reaches the end of the screen
				bullet[i].live = false;  // then it is no longer alive. 
		}
	}
}
void UpdateBullet2(Bullet bullet[], int size)
{
	for (int i = 0; i < size; i++)
	{
		// if the bullet [i] is alive, 
		if(bullet[i].live)
		{
			bullet[i].x -= bullet[i].speed; // continue moving it to the left. note that bullet speed is positive.  
			if(bullet[i].x < 0) // until it reaches the end of the screen
				bullet[i].live = false;  // then it is no longer alive. 
		}
	}
}

// to determine if the bullet from ship collides with ship2 
void CollideBullet(Bullet bullet[], int bSize, SpaceShip &ship2, SpaceShip &ship)
{
	for ( int i = 0; i < bSize; i++)
	{
		if(bullet[i].live)
		{
			if(ship2.lives > 0)
			{		// if bullet is between x and boundx of ship2
					// and between y and boundy of ship2
					// there is a collision. 
				if( bullet[i].x > (ship2.x - ship2.boundx) && 
					bullet[i].x < (ship2.x + ship2.boundx) &&
					bullet[i].y > (ship2.y - ship2.boundy) && 
					bullet[i].y < (ship2.y + ship2.boundy)) 
				{ 
					ship.score++; // increase ship score
					bullet[i].live = false; 
					ship2.lives--; // decrease ship live
				}
			}
		}
	}
}
// to determine if the bullet from ship2 collides with ship
void CollideBullet2(Bullet bullet[], int b2Size, SpaceShip &ship, SpaceShip &ship2)
{
	for ( int i = 0; i < b2Size; i++)
	{
		if(bullet[i].live)
		{
			if(ship.lives > 0)
			{		// if bullet is between x and boundx of ship
					// and btween y and boundy of comet
					// there is a collision. 
				if( bullet[i].x > (ship.x - ship.boundx) && 
					bullet[i].x < (ship.x + ship.boundx) &&
					bullet[i].y > (ship.y - ship.boundy) && 
					bullet[i].y < (ship.y + ship.boundy)) 
				{ 
					ship2.score++; 
					bullet[i].live = false; 
					ship.lives--;
				}
			}
			
		}
	}
}
// this might help with debugging CollideBullet2
void CollideComet(Comet comets[], int cSize, SpaceShip &ship)
{
	for ( int i = 0; i < cSize ; i++)
	{	
		if(comets[i].live)
		{
			if((comets[i].x - comets[i].boundx) < (ship.x + ship.boundx) &&
				(comets[i].x + comets[i].boundx) > (ship.x - ship.boundx) &&
				(comets[i].y - comets[i].boundy) < (ship.y + ship.boundy) &&
				(comets[i].y + comets[i].boundy) > (ship.y - ship.boundy))
			{
				ship.lives--;
				comets[i].live = false;  
			}
		}
	}
}

// Congrats!! you can now make any 2D game with good logic. 
//
//

// */