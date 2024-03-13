## Kurt Hahn
## 08 December 2023
## PHYS 225 Introduction to Computational Physics and Programming
## Final Project - Gravity in 2D


import pygame
import math
import numpy as np


TITLE = "Gravity in 2D"

WINDOWWIDTH, WINDOWHEIGHT = 1200, 600
TIMEDELAY = 1 #milliseconds
DELTAT = 1e2

ACCELERATIONFACTOR = 1
COLLISIONS = False
BOUNCE = 1

EXPONENT = 2
RESETPATHS = True

SHOWCENTEROFMASS = True

WHITE = (255, 255, 255)
RED = (255, 0, 0)
ORANGE = (255, 127, 0)
YELLOW = (240, 240, 0)
GREEN = (0, 255, 0)
BLUE = (0, 100, 255)
PURPLE = (75, 0, 130)
BLACK = (0, 0, 0)
MAGENTA = (180, 29, 245)
GRAY = (128, 128, 128)

CENTER = (WINDOWWIDTH // 2, WINDOWHEIGHT // 2)
YMARGIN = 50
XMARGIN = 50
OVERLAPMARGIN = 9
WINDOWSHIFT = 50
WINDOWZOOM = 1e-8

G = 6.674e-11*ACCELERATIONFACTOR #N*m^2/kg^2

SUNMASS = 2e30 #kg
SUNRADIUS = 7e8 #m

EARTHMASS = 6e24 #kg
EARTHRADIUS = 6.4e6 #m
EARTHDISTANCE = 1.5e11 #m
EARTHORBITALSPEED = 3e4 #m/s

MOONMASS = 7.3e22 #kg
MOONRADIUS = 1.7e6 #m
MOONDISTANCE = 3.8e8 #m
MOONORBITALSPEED = 1e3 #m/s

ISSMASS = 4.5e5 #kg
ISSLENGTH = 91 #m
ISSDISTANCE = 420e3 #m
ISSORBITALSPEED = 7.7e3 #m/s

MARSMASS = 6.42e23
MARSRADIUS = 3.4e6 #m
MARSDISTANCE = 2.28e11 #m
MARSORBITALSPEED = 2.41e4 #m/s

MUSIC = "cornfield_chase.mp3"


class Ball():

    def __init__(self, window, xPos, yPos, xVel, yVel, MASS, RADIUS, COLOR, name):

        self.WINDOW= window
        self.COLOR = COLOR

        self.xPos = xPos
        self.yPos = yPos

        self.xVel = xVel
        self.yVel = yVel

        self.MASS = MASS
        self.RADIUS = RADIUS

        self.name = name


class Simulation():

    def __init__(self, ballset, WINDOW):

        self.ballset = ballset
        self.WINDOW = WINDOW
        self.xShift = 0
        self.yShift = 0
        self.zoom = 1
        self.windowShift = WINDOWSHIFT

        self.indexFocus = 0
        self.deltaT = DELTAT
        self.timeElapsed = 0 #milliseconds

    def return_distance(self, ball1, ball2):
        return math.sqrt((ball1.xPos - ball2.xPos)**2 + (ball1.yPos - ball2.yPos)**2)

    def return_acceleration(self, ball1):

        xAcc = 0
        yAcc = 0

        for ball2 in self.ballset:

            if ball1 == ball2:
                continue

            radius = self.return_distance(ball1, ball2)

            angle = math.atan2(ball2.yPos - ball1.yPos, ball2.xPos - ball1.xPos)

            acc = G * ball2.MASS / (radius**EXPONENT)
            xAcc += acc * math.cos(angle)
            yAcc += acc * math.sin(angle)

        return xAcc, yAcc

    def set_next_frame_velocities(self, ball):

        xAcc, yAcc = self.return_acceleration(ball)

        ball.xVel += xAcc * self.deltaT
        ball.yVel += yAcc * self.deltaT

    def set_next_frame_positions(self, ball):

        ball.xPos += ball.xVel * self.deltaT
        ball.yPos += ball.yVel * self.deltaT

    def check_for_collision(self, ball1, ball2):
        return self.return_distance(ball1, ball2) < (ball1.RADIUS + ball2.RADIUS)
    
    def check_for_complete_overlap(self, ball1, ball2):
        return self.return_distance(ball1, ball2) < OVERLAPMARGIN

    ## This crap came from Wikipedia
    ## https://en.wikipedia.org/wiki/Elastic_collision 
    def set_velocities_after_collision(self, ball1, ball2):

        m1 = ball1.MASS
        m2 = ball2.MASS

        v1 = math.sqrt(ball1.xVel**2 + ball1.yVel**2)
        v2 = math.sqrt(ball2.xVel**2 + ball2.yVel**2)

        # Movement Angles
        theta1 = math.atan2(ball1.yVel, ball1.xVel)
        theta2 = math.atan2(ball2.yVel, ball2.xVel)

        # Contact Angle
        phi = math.atan2(ball2.yPos - ball1.yPos, ball2.xPos - ball1.xPos)

        v1LargeCalc = (v1*math.cos(theta1-phi)*(m1-m2)+2*m2*v2*math.cos(theta2-phi))/(m1+m2)
        v1SmallCalc = v1*math.sin(theta1-phi)

        new_xVel1 = v1LargeCalc*math.cos(phi) + v1SmallCalc*math.cos(phi+np.pi/2)
        new_yVel1 = v1LargeCalc*math.sin(phi) + v1SmallCalc*math.sin(phi+np.pi/2)

        v2LargeCalc = (v2*math.cos(theta2-phi)*(m2-m1)+2*m1*v1*math.cos(theta1-phi))/(m2+m1)
        v2SmallCalc = v2*math.sin(theta2-phi)

        new_xVel2 = v2LargeCalc * math.cos(phi) + v2SmallCalc * math.cos(phi + np.pi / 2)
        new_yVel2 = v2LargeCalc * math.sin(phi) + v2SmallCalc * math.sin(phi + np.pi / 2)

        ball1.xVel = new_xVel1 * BOUNCE
        ball1.yVel = new_yVel1 * BOUNCE
        ball2.xVel = new_xVel2 * BOUNCE
        ball2.yVel = new_yVel2 * BOUNCE

    def set_positions_after_collision(self, ball1, ball2):
        
        distance = self.return_distance(ball1, ball2)

        overlap = (ball1.RADIUS + ball2.RADIUS) - distance

        angle = math.atan2(ball2.yPos - ball1.yPos, ball2.xPos - ball1.xPos)

        ball1.xPos -= (overlap / 2) * math.cos(angle)
        ball1.yPos -= (overlap / 2) * math.sin(angle)

        ball2.xPos += (overlap / 2) * math.cos(angle)
        ball2.yPos += (overlap / 2) * math.sin(angle)

    def draw_balls(self):

        for ball in self.ballset:

            imageX, imageY = self.calculate_image_position(ball.xPos, ball.yPos)
            imageRadius = ball.RADIUS * self.zoom

            pygame.draw.circle(self.WINDOW, ball.COLOR, (imageX, imageY), imageRadius)

    def evolve(self):
             
        if not COLLISIONS:

            for ball in self.ballset:
                self.set_next_frame_velocities(ball)

            for ball in self.ballset:
                self.set_next_frame_positions(ball)
        
        else: 
        
            ## Collision mechanics
            for index1, ball1, in enumerate(self.ballset):

                isStuck = False

                for index2, ball2 in enumerate(self.ballset):

                    if index1 <= index2:
                        continue

                    if not self.check_for_collision(ball1, ball2):
                        continue
                    
                    if self.check_for_complete_overlap(ball1, ball2):
                        isStuck = True
                        continue
                   

                    self.set_velocities_after_collision(ball1, ball2)
                    
                ## Gravity Mechanics
                if isStuck:
                    ball1.xVel = 0
                    ball2.yVel = 0
                else:
                    self.set_next_frame_velocities(ball1)

            for ball in self.ballset:
                self.set_next_frame_positions(ball)

            ## Rectify overlaps
            for index2, ball2 in enumerate(self.ballset):
                if index1 <= index2:
                    continue

                if not self.check_for_collision(ball1, ball2):
                    continue

                self.set_positions_after_collision(ball1, ball2)

    def return_center_of_mass(self):

        xResult = 0
        yResult = 0
        totalMass = 0

        for ball in self.ballset:

            xResult += ball.MASS * ball.xPos
            yResult += ball.MASS * ball.yPos
            totalMass += ball.MASS

        xResult /= totalMass
        yResult /= totalMass

        return xResult, yResult
    
    def draw_center_of_mass_and_lines(self):

        xPos, yPos = self.return_center_of_mass()

        imageX, imageY = self.calculate_image_position(xPos, yPos)
        imageRadius = self.zoom
               
        pygame.draw.circle(self.WINDOW, WHITE, (imageX, imageY), imageRadius)

        for ball in self.ballset:

            imageX2, imageY2 = self.calculate_image_position(ball.xPos, ball.yPos)

            pygame.draw.line(self.WINDOW, WHITE, (imageX, imageY), (imageX2, imageY2))

    def calculate_image_position(self, xPos, yPos):

        focusBall = self.ballset[self.indexFocus]

        imageX = CENTER[0] + (xPos - self.xShift - focusBall.xPos) * self.zoom
        imageY = CENTER[1] - (yPos - self.yShift - focusBall.yPos) * self.zoom 

        return imageX, imageY
 

class Main():

    def main():

        pygame.init()

        WINDOW = pygame.display.set_mode((WINDOWWIDTH, WINDOWHEIGHT))
        pygame.display.set_caption(TITLE)
        focusTextFont = pygame.font.Font("freesansbold.ttf", 32)
        labelTextFont = pygame.font.Font("freesansbold.ttf", 8)
        deltatTextFont = pygame.font.Font("freesansbold.ttf", 20)
        commandsTextFont = pygame.font.Font("freesansbold.ttf", 15)
        timeElapsedTextFont = pygame.font.Font("freesansbold.ttf", 20)
        icon = pygame.image.load('image1.png')
        pygame.display.set_icon(icon)

        pygame.mixer.music.load(MUSIC)
        pygame.mixer.music.set_volume(1)
        pygame.mixer.music.play(-1)

 
        Sun = Ball(WINDOW, 0, 0, 0, 0, SUNMASS, SUNRADIUS, YELLOW, "Sun")
        Earth = Ball(WINDOW, 0, EARTHDISTANCE, EARTHORBITALSPEED, 0, EARTHMASS, EARTHRADIUS, BLUE, "Earth")
        Moon = Ball(WINDOW, 0, EARTHDISTANCE+MOONDISTANCE, EARTHORBITALSPEED+MOONORBITALSPEED, 0, MOONMASS, MOONRADIUS, WHITE, "Moon")
        InternationalSpaceStation = Ball(WINDOW, 0, EARTHDISTANCE+EARTHRADIUS+ISSDISTANCE, EARTHORBITALSPEED+ISSORBITALSPEED, 0, ISSMASS, ISSLENGTH, GRAY, "ISS") 
        Mars = Ball(WINDOW, 0, MARSDISTANCE, MARSORBITALSPEED, 0, MARSMASS, MARSRADIUS, RED, "Mars")
        ballset = [Sun, Earth, Moon, InternationalSpaceStation, Mars]


        simulation = Simulation(ballset, WINDOW)
        simulation.zoom = 1e-6


        run = True
        pause = True
        inputTick = 0
        while run:

            pygame.time.delay(TIMEDELAY)


            ## Pressed input actions
            pressedKeys = pygame.key.get_pressed()
            if inputTick == 0:
                if pressedKeys[pygame.K_LSHIFT] or pressedKeys[pygame.K_RSHIFT]:
                    WINDOWZOOM = 3e-8
                if pressedKeys[pygame.K_UP]:
                    WINDOW.fill(BLACK)
                    simulation.yShift += simulation.windowShift / simulation.zoom
                if pressedKeys[pygame.K_DOWN]:
                    WINDOW.fill(BLACK)
                    simulation.yShift -= simulation.windowShift / simulation.zoom
                if pressedKeys[pygame.K_RIGHT]:
                    WINDOW.fill(BLACK)
                    simulation.xShift += simulation.windowShift / simulation.zoom
                if pressedKeys[pygame.K_LEFT]:
                    WINDOW.fill(BLACK)
                    simulation.xShift -= simulation.windowShift / simulation.zoom
                if pressedKeys[pygame.K_p]:
                    WINDOW.fill(BLACK)
                    simulation.zoom += WINDOWZOOM
                if pressedKeys[pygame.K_m] and simulation.zoom - WINDOWZOOM >= 0:
                    WINDOW.fill(BLACK)
                    simulation.zoom -= WINDOWZOOM
                
                inputTick = 10
                WINDOWZOOM = 1e-8
            else:
                inputTick -= 1


            ## Clicked input actions
            for event in pygame.event.get():
                if event.type == pygame.QUIT:
                    run = False
                if event.type == pygame.KEYDOWN:
                    if event.key == pygame.K_SPACE:
                        if pause == False: pause = True
                        else: pause = False
                    if event.key == pygame.K_ESCAPE:
                        pygame.quit()
                    if event.key == pygame.K_LEFTBRACKET:
                        simulation.indexFocus -= 1
                        simulation.indexFocus %= len(simulation.ballset)
                        simulation.xShift = 0
                        simulation.yShift = 0
                        WINDOW.fill(BLACK)
                    if event.key == pygame.K_RIGHTBRACKET:
                        simulation.indexFocus += 1
                        simulation.indexFocus %= len(simulation.ballset)
                        simulation.xShift = 0
                        simulation.yShift = 0
                        WINDOW.fill(BLACK)
                    if event.key == pygame.K_1:
                        simulation.deltaT = 1e-3
                    if event.key == pygame.K_2:
                        simulation.deltaT = 1e-2
                    if event.key == pygame.K_3:
                        simulation.deltaT = 1e-1
                    if event.key == pygame.K_4:
                        simulation.deltaT = 1
                    if event.key == pygame.K_5:
                        simulation.deltaT = 1e1
                    if event.key == pygame.K_6:
                        simulation.deltaT = 1e2
                    if event.key == pygame.K_7:
                        simulation.deltaT = 1e3
                    if event.key == pygame.K_8:
                        simulation.deltaT = 1e4
                    if event.key == pygame.K_9:
                        simulation.deltaT = 1e5
                    if event.key == pygame.K_0:
                        pause = True


            ## Reset Screen
            if RESETPATHS: WINDOW.fill(BLACK)


            ## Draw center of mass and lines
            if SHOWCENTEROFMASS: simulation.draw_center_of_mass_and_lines()


            ## Draw balls
            simulation.draw_balls()


            ## Add corner text
            focusText = focusTextFont.render(f"Focusing on: {simulation.ballset[simulation.indexFocus].name}", True, GREEN, BLUE)
            WINDOW.blit(focusText, (0,0))


            ## Add labels for each body
            for ball in simulation.ballset:
                
                imageX, imageY = simulation.calculate_image_position(ball.xPos, ball.yPos)

                text = labelTextFont.render(f"{ball.name}", True, WHITE, None)
                WINDOW.blit(text, (imageX, imageY))


            ## Add deltaT text box
            deltatText = deltatTextFont.render(f"Δt = {simulation.deltaT} ({simulation.deltaT * 1000}x)", True, GREEN, BLUE)
            WINDOW.blit(deltatText, (0, 32))


            ## Add text box for key commands
            scrollCommandText = commandsTextFont.render("Use arrow keys to scroll", True, WHITE, None)
            zoomCommandText = commandsTextFont.render("Press \"p\" to zoom in and \"m\" to zoom out", True, WHITE, None)
            pauseCommandText = commandsTextFont.render("Press SPACE to pause", True, WHITE, None)
            fastForwardCommandText = commandsTextFont.render("Use number keys to change Δt (1-6 recommended)", True, WHITE, None)
            focusCommandText = commandsTextFont.render("Press \"[\" or \"]\" to change reference frame", True, WHITE, None)
            terminateCommandText = commandsTextFont.render("Press ESC to terminate simulation", True, WHITE, None)
            commandList = [scrollCommandText, zoomCommandText, pauseCommandText, fastForwardCommandText, focusCommandText, terminateCommandText]
            for index, command in enumerate(commandList):
                WINDOW.blit(command, (0, 100+20*index))


            ## Add text box for Time Elapsed
            timeElapsedMillisecondsText = timeElapsedTextFont.render(f"t = {int(simulation.timeElapsed)} ms", True, GREEN, BLUE)
            WINDOW.blit(timeElapsedMillisecondsText, (0, 52))
            timeElapsedYearsText = timeElapsedTextFont.render(f"t = {int(simulation.timeElapsed / 31558149763.5456)} years", True, GREEN, BLUE)
            WINDOW.blit(timeElapsedYearsText, (0, 72))


            ## Update screen
            pygame.display.update()
            

            ## Skip calculate next frame if paused
            if pause: continue


            ## Calculate next frame
            simulation.evolve()
            simulation.timeElapsed += simulation.deltaT * 1000


        pygame.quit()


if __name__ == "__main__":
    Main.main()