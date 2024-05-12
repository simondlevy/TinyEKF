SKETCH = SensorFusion

PORT = /dev/ttyACM0

BOARD = arduino:avr:uno

LIBS = ../../..

build: $(SKETCH).ino 
	arduino-cli compile  --libraries $(LIBS) --fqbn $(BOARD) $(SKETCH).ino

flash:
	arduino-cli upload -p $(PORT) --fqbn $(BOARD) .

listen:
	miniterm $(PORT) 115200

edit:
	vim $(SKETCH).ino
