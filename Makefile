.PHONY=all clean bonsai update

all:
	+$(MAKE) -C bonsai

clean:
	+$(MAKE) clean -C bonsai

unit:
	+$(MAKE) unit -C bonsai

bonsai:
	+$(MAKE) bonsai -C bonsai

update:
	+$(MAKE) update -C bonsai
