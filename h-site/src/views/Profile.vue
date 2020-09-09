<template>
	<div>
		<div v-if="auth">
			<p>Email: {{ user.email }}</p>
			<p v-if="user.displayName">Username: {{ user.displayName }}</p>
			<div v-if="!user.displayName">Your username is not set:
				<p>
					<input v-model="username" placeholder="Update not available at the moment">
					<button id="setname" @click="setUsername">Set username</button></p>
			</div>
			<button v-if="!verify" id="delete" @click="verifyDelete">Delete account</button>
			<div v-if="verify">
				Are you sure you want to delete your account? 
				<p>
					<button id="delete" @click="deleteAccount">Yes</button>
					<button id="delete" @click="verifyDelete">No, I changed my mind</button>
				</p>
			</div>
		</div>
		<div v-if="!auth">You need to sign in</div>
	</div>
</template>

<script>

import { mapGetters, mapActions } from 'vuex'

export default {
	name: "Profile",
	data() {
		return {
			username: '',
			verify: false
		}
	},
	methods: {
		...mapActions(['updateDisplayName', 'deleteAccount']),
		setUsername() {
			this.updateDisplayName(this.username)
			console.log("usernameset: ", this.username)
		},
		verifyDelete() {
			this.verify = !this.verify
		}
	},
	computed: {
		...mapGetters(['auth', 'user'])
	}
}
</script>

<style>
#setname {
	border-style: outset;
	background-color: silver;
	border-radius: 5px;
}
#delete {
	border-style: inset;
	background-color: red;
	border-radius: 5px;
}
</style>