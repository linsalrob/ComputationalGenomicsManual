# Upgrading Images

We strive to keep the SDSU Computational Genomics Image upto date, and sometimes that means updating the system software, and sometimes fixing installation bugs and issues. From time to time you may need to upgrade the image.

There are a couple of steps that can make life easier for you!.

If you have a running image with some data on it, leave it running. Now, launch a new instance of the latest version of our image. You now have two running instances, an old instance with your data and a shiny new instance with no data on it.

We're going to copy the data across, but before we do that we need to move the private [ssh key](Linux#public-and-private-ssh-keys) onto the server (but don't worry, because we're going to delete the server shortly).

scp your `pem` file from your laptop to the server.

```bash
scp ~/.ssh/id_rsa.pem ec2-user@xxx.xxx.xxx.xxx:~/.ssh/
```
_Note_ You need to change ~/.ssh/id_rsa.pem to the location of your `pem` file that you acknowledge that you have when you start an AWS instance. Also, you need to change `xxx.xxx.xxx.xxx` to the IP of the older server that still has your data on it.

Next, login to the older server

```bash
ssh xxx.xxx.xxx.xxx
```

and then scp the data to the new server

```bash
scp -r * ec2-user@yyy.yyy.yyy.yyy:
```
_Note_ change `yyy.yyy.yyy.yyy` to the ip address of the new server. Also, don't forget the colon on the end of the line!

This will ask you if you want to accept the identity of the new server. Go ahead and type `yes`.

If you login to the new server, you should see all your data there! If not, take a look at the old server and see what you did wrong ... did you forget the `:` on the end of the line?

Once you have all the data on the new server, you can log back into the AWS console. Select the older server and choose `Actions --> Instance State --> Terminate`. That will delete it, and everything including your private key.
