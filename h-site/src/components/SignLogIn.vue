<template>
	<div id="mark">
		<div v-if="(!login || !signup) && !auth">
			<button class="login" @click="setlogin">
				Logg inn
			</button>
			<button class="signup" @click="setsignup">
				Sign up
			</button>
		</div>
		<div v-if="auth">
			<p v-if="user.displayName">{{ user.displayName }}</p>
			<p v-else>{{ user.email }}</p>
			<button class="login" @click="signOut">Sign out</button>
		</div>

		<div id="login-box" v-if="login">
			Log inn with your email
			<p>Email: <input v-model="usremailL" placeholder="email@address.com"></p>
			
			<p>Password: <input @keyup.enter="submitLogin" type="password" v-model="passwordL" placeholder="password"></p>

			<p>
				<button id="con" @click="submitLogin">Submit</button>
				<button id="con" @click="login = false">Cancel</button>
			</p>
		</div>

		<div id="login-box" v-if="signup">
			Sign up with your email
			<p>Email: <input v-model="usremailS" placeholder="email@address.com"></p>
			
			<p>Password: <input type="password" v-model="passwordS" placeholder="password"></p>

			<p>
				<button id="con" @click="submitSignup">Submit</button>
				<button id="con" @click="signup = false">Cancel</button>
			</p>
		</div>
	</div>
</template>

<script>
import { mapActions, mapGetters } from 'vuex'

	export default {
		name: 'SignLogIn',
		data() {
			return {
				login: false,
				signup: false,
				passwordL: '',
				usremailL: '',
				passwordS: '',
				usremailS: ''
			}
		},
		methods: {
			...mapActions(['authentication', 'logOut', 'signUp']),
			setlogin() {
				this.login = !this.login
				this.signup = false
			},
			submitLogin() {
				//console.log(this.passwordL, this.usremailL)
				this.authentication({
					email: this.usremailL, 
					password: this.passwordL})
				this.login = !this.login
			},
			setsignup() {
				this.signup = !this.signup
				this.login = false
			},
			submitSignup() {
				//console.log("signup")
				this.signUp({
					email: this.usremailS,
					password: this.passwordS
				})
				this.signup = !this.signup
			},
			signOut() {
				this.logOut()
			}
		},
		computed: {
			...mapGetters(['auth', 'user'])
		}
	}
</script>

<style> 
	button {
		border: none;
		background-color: black;
	}
	.login {
		color: red;
	}
	.signup {
		color: blue;
	}
	#mark {
		position: relative;
		right: -80%;
		width: 20%;
		top: -25px;
	}
	#login-box {
		position: fixed;
		right: 0%;
		background: snow;
		border-style: outset;
		padding: 5px;
		text-align: center;
		width: 50%
	}
	#con {
		border-radius: 5px;
		border-style: inset;
		background-color: lightblue;

	}
</style>
